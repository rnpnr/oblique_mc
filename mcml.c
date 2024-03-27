/* Oblique Monte Carlo
 * Written by Randy Palamar - March 2024
 * Based on: Hybrid-MC by Li Li
 * Refrence: Wang, "Biomedical optics", chpater 3 & 6
 */

/* Details:
 * - Uses Cartesian coordinates.
 * - Beam is incident in x-z plane.
 */

#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define SGN(x) ((x) >= 0 ? 1 : -1)
#define ABS(x) ((x) >= 0 ? x : -x)

#define IS_SPACE(x) ((x) == '\n' || (x) == '\r' || (x) == '\t' || (x) == ' ')

#define ZERO   1.0e-6
#define BIGVAL 1.0e6
#define EPS    1.0e-9

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef uint8_t   u8;
typedef double    f64;
typedef ptrdiff_t size;

typedef struct { u8 *data; size len; } s8;

typedef struct {
	f64 x, y ,z;
} Vec3;

typedef struct {
	u32 Nx, Ny;
	f64 *b;
} Mat2;

typedef struct {
	Vec3 pos;
	Vec3 dir;
	f64 w;
	u32 n_scatters;
	u32 dead;
} Photon;

/* global data that will not be modified after startup */
static struct {
	f64 dx, dy;
	f64 xoff, yoff;
	u32 Nx, Ny;
	u32 N_photons;

	f64 mu_a, mu_s, mu_t;
	f64 g, d;

	f64 n, n0;
	f64 theta_i;
} gctx;

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

static s8
read_file(const char *fname)
{
	int fd;
	s8 ret;

	if ((fd = open(fname, O_RDONLY)) < 0)
		die("open: %s\n", fname);

	ret.len = lseek(fd, 0L, SEEK_END);
	ret.data = malloc(ret.len);
	if (ret.data == NULL)
		die("malloc\n");
	lseek(fd, 0L, SEEK_SET);

	size rlen;
	do {
		rlen = read(fd, ret.data, ret.len);
	} while (rlen == -1 && errno == EINTR);
	close(fd);

	if (rlen != ret.len)
		die("read: %s\n", fname);

	return ret;
}

static s8
split_at_space(s8 s)
{
	s8 ret;

	/* skip leading white space */
	while (IS_SPACE(*s.data) && s.len) {
		s.data++;
		s.len--;
	}
	ret.data = s.data;

	while (!IS_SPACE(*s.data) && s.len) {
		s.data++;
		s.len--;
	}
	s.data[0] = 0;
	ret.len = s.data - ret.data;
	return ret;
}

/* fills in global ctx based on input */
static void
load_input(const char *name)
{
	s8 f = read_file(name);
	s8 s = split_at_space(f);

	char *e;
	#define ADVANCE_S(_s) do {\
		_s.data = (u8 *)e + 1;\
		_s.len = f.len - (f.data - s.data);\
		_s = split_at_space(_s);\
	} while (0);

	errno = 0;
	gctx.N_photons = (u64)strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.theta_i = strtod((char *)s.data, &e) * M_PI / 180;
	ADVANCE_S(s);
	gctx.dx = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.dy = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.Nx = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.Ny = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.n = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.mu_a = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.mu_s = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.g = strtod((char *)s.data, &e);
	ADVANCE_S(s);
	gctx.d = strtod((char *)s.data, &e);

	gctx.xoff = gctx.dx * gctx.Nx / 2;
	gctx.yoff = gctx.dy * gctx.Ny / 2;
	gctx.mu_t = gctx.mu_a + gctx.mu_s;

	gctx.n0 = 1.0;

	#undef ADVANCE_S

	if (errno != 0)
		die("strtod: errno = %d\n", errno);

	free(f.data);
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

static void
launch_photon(Photon *p)
{
	p->pos = (Vec3){0, 0, 0};
	p->dir.x = sin(gctx.theta_i) * gctx.n0 / gctx.n;
	p->dir.y = 0;
	p->dir.z = sqrt(1 - p->dir.x * p->dir.x);
	p->n_scatters = 0;
	p->dead = 0;

	f64 sin_ti = sin(gctx.theta_i);
	f64 cos_ti = cos(gctx.theta_i);
	f64 sin_tt = p->dir.x;
	f64 cos_tt = p->dir.z;
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
				u32 ri = (p->pos.x + gctx.xoff) / gctx.dx;
				u32 ci = (p->pos.y + gctx.yoff) / gctx.dy;
				Rd_xy.b[ri * Rd_xy.Ny + ci] += p->w;
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

	launch_photon(p);
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
	if (argc != 3)
		die("usage: %s input_file output_file\n", argv[0]);

	load_input(argv[1]);
	char *outfile = argv[2];

	random_init();

	alloc_mat2(&Rd_xy, gctx.Nx, gctx.Ny);

	/* Propagate Photons */
	time_t tstart, tend;
	time(&tstart);
	Photon p;
	for (u32 i = 1; i <= gctx.N_photons; i++) {
		simulate_photon(&p);
		if (i % (gctx.N_photons / 10) == 0)
			printf("[%u/%u] photons done!\n", i, gctx.N_photons);
	}
	time(&tend);
	printf("Simulation took: %ld [s]\n", tend - tstart);


	FILE *fd = fopen(outfile, "w");
	if (fd == NULL)
		die("rip can't open output file: %s\n", outfile);

	fputs("/////////////Grid coordinates///////////\nx (cm): ", fd);
	for (u32 i = 0; i < gctx.Nx; i++)
		fprintf(fd, "%e ", (i + 0.5) * gctx.dx - gctx.xoff);
	fputs("\ny (cm): ", fd);
	for (u32 i = 0; i < gctx.Nx; i++)
		fprintf(fd, "%e ", (i + 0.5) * gctx.dy - gctx.yoff);
	fputs("\n/////////////////////////////////////////\n", fd);
	fputs("///////////////Rd_xy/////////////////////\n", fd);

	f64 scale = gctx.N_photons * gctx.dx * gctx.dy;
	f64 *b = Rd_xy.b;
	for (u32 i = 0; i < Rd_xy.Nx; i++) {
		for (u32 j = 0; j < Rd_xy.Ny; j++)
			fprintf(fd, "%e ", b[j] / scale);
		fputc('\n', fd);
		b += Rd_xy.Ny;
	}
	fputs("\n/////////////////////////////////////////\n", fd);
	fclose(fd);

	return 0;
}
