/* Oblique Monte Carlo
 * Written by Randy Palamar - March-April 2024
 * Based on: Hybrid-MC by Li Li
 * Reference: Wang, "Biomedical optics", chapter 3 & 6
 */

/* Details:
 * - Mostly uses Cartesian coordinates. Polar coordinates are used for
 *   specifying incident location.
 * - Beam faces in positive z direction incident from anywhere in r-theta
 *   plane. Initial launch direction is always towards origin.
 */

#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SGN(x)     ((x) >= 0 ? 1 : -1)
#define ABS(x)     ((x) >= 0 ? x : -x)
#define LEN(a)     (sizeof(a) / sizeof(*a))
#define DEG2RAD(a) ((a) * M_PI / 180.0)

#define ZERO   1.0e-6
#define BIGVAL 1.0e6
#define EPS    1.0e-9

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef uint8_t   u8;
typedef double    f64;
typedef ptrdiff_t size;

typedef struct { u8 *data; size len; } s8;
#define s8(s) (s8){(u8 *)s, (LEN(s) - 1)}

typedef struct { f64 x, y ,z; }               Vec3;
typedef struct { f64 r, theta, z; }           Vec3Pol;
typedef struct { f64 top, bot, left, right; } Rect;
typedef struct { u32 Nx, Ny; f64 *b; }        Mat2;

typedef struct {
	Vec3 pos;
	Vec3 dir;
	f64 w;
	u32 n_scatters;
	u32 dead;
} Photon;

#include "config.h"

/* these will be modified; FIXME: multithreading */
static Mat2 Rd_xy;
static u64 rand_state[2];

static void
die(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);

	exit(1);
}

static Vec3
polar_to_cartesian(Vec3Pol v)
{
	return (Vec3){
		.x = v.r * cos(v.theta),
		.y = v.r * sin(v.theta),
		.z = v.z
	};
}

static Vec3Pol
cartesian_to_polar(Vec3 v)
{
	return (Vec3Pol){
		.r = sqrt(v.x * v.x + v.y * v.y),
		.theta = atan2(v.y, v.x),
		.z = v.z
	};
}

/* concatenates nstrs together. 0 terminates output for use with bad APIs */
static s8
s8concat(s8 *strs, size nstrs)
{
	s8 out = {0};
	for (size i = 0; i < nstrs; i++)
		out.len += strs[i].len;
	out.data = malloc(out.len + 1);
	if (out.data == NULL)
		die("s8concat\n");
	for (size i = 0, off = 0; i < nstrs; i++) {
		if (!memcpy(out.data + off, strs[i].data, strs[i].len))
			die("memcpy\n");
		off += strs[i].len;
	}
	out.data[out.len] = 0;
	return out;
}

static void
dump_output(s8 pre)
{
	s8 xy = s8("_xy.tsv");
	s8 rd = s8("_Rd_xy.csv");
	s8 cat[2] = { pre, xy };
	s8 out = s8concat(cat, 2);

	FILE *fd = fopen((char *)out.data, "w");
	if (fd == NULL)
		die("rip can't open output file: %s\n", out.data);
	fputs("x [cm]\ty [cm]\n", fd);
	for (u32 i = 0; i < gctx.Nx; i++) {
		f64 x = (i + 0.5) * gctx.dx - gctx.xoff;
		f64 y = (i + 0.5) * gctx.dy - gctx.yoff;
		fprintf(fd, "%e\t%e\n", x, y);
	}
	fclose(fd);

	free(out.data);
	cat[1] = rd;
	out = s8concat(cat, 2);

	fd = fopen((char *)out.data, "w");
	if (fd == NULL)
		die("rip can't open output file: %s\n", out.data);

	f64 scale = gctx.N_photons * gctx.dx * gctx.dy;
	f64 *b = Rd_xy.b;
	for (u32 i = 0; i < Rd_xy.Nx; i++) {
		for (u32 j = 0; j < Rd_xy.Ny; j++)
			fprintf(fd, "%e,", b[j] / scale);
		fseek(fd, -1, SEEK_CUR);
		fputc('\n', fd);
		b += Rd_xy.Ny;
	}
	fclose(fd);
	free(out.data);
}

/* fill in extra global ctx values */
static void
init(void)
{
	if (gctx.Nx != gctx.Ny)
		die("Nx != Ny, output must be square!\n");
	f64 w = gctx.extent.right - gctx.extent.left;
	f64 h = gctx.extent.top - gctx.extent.bot;
	gctx.dx = w / gctx.Nx;
	gctx.dy = h / gctx.Ny;
	gctx.xoff = w / 2;
	gctx.yoff = h / 2;
	gctx.mu_t = gctx.mu_a + gctx.mu_s;
}

static void
random_init(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	rand_state[0] = (intptr_t)&printf ^ ts.tv_sec
	              ^ ((u64)ts.tv_nsec * 0xAC5533CD);
	rand_state[1] = (intptr_t)&malloc ^ ts.tv_sec
	              ^ ((u64)ts.tv_nsec * 0xAC5533CD);
}

static u64
xoroshiro128plus(u64 s[2])
{
	u64 s0 = s[0];
	u64 s1 = s[1];
	u64 result = s0 + s1;
	s1 ^= s0;
	s[0] = ((s0 << 24) | (s0 >> 40)) ^ s1 ^ (s1 << 16);
	s[1] = (s1 << 37) | (s1 >> 27);
	return result;
}

static f64
random_uniform(void)
{
	return xoroshiro128plus(rand_state) / (f64)UINT64_MAX;
}

static void
alloc_mat2(Mat2 *m, u32 x, u32 y)
{
	m->b = calloc(x * y, sizeof(f64));
	if (m->b == NULL)
		die("calloc\n");
	m->Nx = x;
	m->Ny = y;
}

static Vec3
launch_direction_from_polar_angle(f64 theta)
{
	f64 rdir = sin(gctx.theta_i) * gctx.n0 / gctx.n;
	return (Vec3){
		.x = -rdir * cos(theta),
		.y = -rdir * sin(theta),
		.z = sqrt(1 - rdir * rdir)
	};
}

static void
launch_photon(Photon *p, Vec3Pol pos)
{
	p->pos = polar_to_cartesian(pos);
	p->dir = launch_direction_from_polar_angle(pos.theta);
	p->n_scatters = 0;
	p->dead = 0;

	Vec3Pol dir = cartesian_to_polar(p->dir);
	f64 sin_ti = sin(gctx.theta_i);
	f64 cos_ti = cos(gctx.theta_i);
	f64 sin_tt = dir.r;
	f64 cos_tt = dir.z;
	f64 raux1 = sin_ti * cos_tt - cos_ti * sin_tt;
	raux1 *= raux1;
	f64 raux2 = sin_ti * cos_tt + cos_ti * sin_tt;
	raux2 *= raux2;
	f64 rsp = raux1 / raux2 + raux1 * (1 - raux2) / raux2 / (1 - raux1);
	rsp *= 0.5;
	p->w = (1.0 - rsp);
}

static void
move_photon(Photon *p, f64 s)
{
	p->pos.x += s * p->dir.x;
	p->pos.y += s * p->dir.y;
	p->pos.z += s * p->dir.z;
}

static void
absorb_photon(Photon *p)
{
	f64 dw = p->w * gctx.mu_a / gctx.mu_t; /* eq 3.24 */
	p->w -= dw;
	if (p->w < 0)
		p->dead = 1;
}

static void
reflect_or_transmit_photon(Photon *p)
{
	f64 sin_ai = sqrt(1 - p->dir.z * p->dir.z);
	f64 sin_at = sin_ai * gctx.n / gctx.n0; /* eq. 3.35 */

	f64 r_ai;
	if (sin_at < ZERO) { /* roughly normal */
		r_ai = (gctx.n - gctx.n0) / (gctx.n + gctx.n0);
		r_ai *= r_ai;
	} else if (sin_at > 1.0) { /* TIR */
		r_ai = 1;
	} else {
		f64 A = sin_ai * sqrt(1 - sin_at * sin_at)
		        - sin_at * sqrt(1 - sin_ai * sin_ai);
		A *= A;
		f64 B = sin_ai * sqrt(1 - sin_at * sin_at)
		        + sin_at * sqrt(1 - sin_ai * sin_ai);
		B *= B;
		r_ai = 0.5 * (2 * A - A * A - A * B) / (B * (1 - A)); /* eq 3.36 */
	}

	f64 r = random_uniform();
	if (r <= r_ai) {
		/* rebound to current layer */
		p->dir.z = -p->dir.z;
	} else {
		/* cross the boundary */
		if (p->dir.z < -ZERO) {
			/* hit upper boundary. record reflactance if
			 * we are in bounds of output region */
			if (p->pos.x > -gctx.xoff && p->pos.x < gctx.xoff &&
			    p->pos.y > -gctx.yoff && p->pos.y < gctx.yoff) {
				u32 ri = (p->pos.y + gctx.yoff) / gctx.dy;
				u32 ci = (p->pos.x + gctx.xoff) / gctx.dx;
				Rd_xy.b[ri * Rd_xy.Nx + ci] += p->w;
			}
		}
		p->dead = 1;
	}
}

static void
scatter_photon(Photon *p)
{
	if (p->dead)
		return;

	f64 cos_t, fei;
	f64 g = gctx.g;
	if (g != 0) {
		f64 r = random_uniform();
		f64 aa = (1 - g * g) / (1 - g + 2 * g * r);
		cos_t = (1 + g * g - aa * aa) / (2 * g);
		if (cos_t < -1)
			cos_t = -1;
		else if (cos_t > 1)
			cos_t = 1;
	} else {
		cos_t = 2 * random_uniform() - 1;
	} /* eq. (3.28) */

	f64 sin_t, sin_fei, cos_fei;
	sin_t = sqrt(1 - cos_t * cos_t);
	fei = 2 * M_PI * random_uniform(); /* eq. (3.29) */
	cos_fei = cos(fei);
	sin_fei = sin(fei);

	f64 uux, uuy, uuz;
	if (ABS(p->dir.z) <= (1 - ZERO)) {
		/* eq. (3.24) */
		Vec3 u = p->dir;
		uux = sin_t * (u.x * u.z * cos_fei - u.y * sin_fei)
		      / sqrt(1 - u.z * u.z) + u.x * cos_t;
		uuy = sin_t * (u.y * u.z * cos_fei + u.x * sin_fei)
		      / sqrt(1 - u.z * u.z) + u.y * cos_t;
		uuz = -sin_t * cos_fei * sqrt(1 - u.z * u.z) + u.z * cos_t;
	} else {
		/* close to normal propagation */
		/* eq. (3.30) */
		uux = sin_t * cos_fei;
		uuy = sin_t * sin_fei;
		uuz = SGN(p->dir.z) * cos_t;
	}
	p->dir.x = uux;
	p->dir.y = uuy;
	p->dir.z = uuz;
	p->n_scatters++;

}

static void
check_photon_life(Photon *p)
{
	if (p->dead)
		return;

	if (p->w > 1e-4)
		return;

	f64 m = 10;
	f64 e = random_uniform();
	if (m * e > 1) {
		p->dead = 1;
	} else {
		p->w *= m;
	}
}

static f64
next_step(f64 s)
{
	if (s < ZERO) {
		f64 r = random_uniform();
		s = -log(r + EPS);
	}
	return s;
}

static f64
step_towards_boundary(Photon *p)
{
	f64 db;
	if (p->dir.z < -ZERO) /* travel up */
		db = -p->pos.z / p->dir.z;
	else if (p->dir.z > ZERO) /* travel down */
		db = (gctx.d - p->pos.z) / p->dir.z;
	else
		db = BIGVAL;
	return db;
}

static void
simulate_photon(Photon *p)
{
	f64 step = 0, boundary_dist = 0;
	do {
		step = next_step(step);
		boundary_dist = step_towards_boundary(p);
		if (boundary_dist * gctx.mu_t <= step) {
			move_photon(p, boundary_dist);
			step -= boundary_dist * gctx.mu_t;
			reflect_or_transmit_photon(p);
		} else {
			move_photon(p, step / gctx.mu_t);
			absorb_photon(p);
			scatter_photon(p);
		}
		check_photon_life(p);
	} while (!p->dead);
}

int
main(int argc, char *argv[])
{
	if (argc != 2)
		die("usage: %s output_prefix\n", argv[0]);
	s8 pre = (s8){.data = (u8 *)argv[1], .len = strlen(argv[1])};

	init();
	random_init();

	alloc_mat2(&Rd_xy, gctx.Nx, gctx.Ny);

	/* Propagate Photons */
	time_t tstart, tend;
	time(&tstart);

	/* cache starting photon; nothing at launch changes between runs */
	Photon p_start;
	launch_photon(&p_start, gctx.incidence_location);
	for (u32 i = 1; i <= gctx.N_photons; i++) {
		/* Photon is 64 bytes. this will use SIMD if available.
		 * otherwise compiler will just insert a memcpy call here */
		Photon p = p_start;
		simulate_photon(&p);
		if (i % (gctx.N_photons / 10) == 0)
			printf("[%u/%u] photons done!\n", i, gctx.N_photons);
	}
	time(&tend);
	printf("Simulation took: %ld [s]\n", tend - tstart);

	dump_output(pre);

	return 0;
}
