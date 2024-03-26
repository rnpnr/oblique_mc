/* Providing global variables and functions*/
#include <math.h>
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;

//global constants
#define PI M_PI

#ifndef ZERO
	#define ZERO 1.0e-6 //for small deflection angle
#endif

#ifndef n0
	#define n0 1.0 //refractive index of ambient medium
#endif

#ifndef BIGVAL
	#define BIGVAL 1e6//big value
#endif

#define EPS 1.0e-9

//global variables
#ifdef GLOBALORIGIN
	#define GLOBAL
#else
	#define GLOBAL extern
#endif
GLOBAL double	dx;
GLOBAL double	dy;
GLOBAL int	Nx;
GLOBAL int	Ny;
GLOBAL double	zc;
GLOBAL double	n;
GLOBAL double	mua;
GLOBAL double	mus;
GLOBAL double	mut;
GLOBAL double	g;
GLOBAL double	d;

GLOBAL double	th_in; // incident angle in degree

GLOBAL double	x_offset;
GLOBAL double	y_offset;

struct Mat {
	u32 Nx, Ny;
	double *b;
};

GLOBAL struct Mat Rd_xy;

#define SIGN(x) ((x) >= 0 ? 1 : -1)

void random_init(void);
double random_uniform(void);

void alloc_mat(struct Mat *, int, int);
