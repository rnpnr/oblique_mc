/* global data that will not be modified at runtime */

/* optional prefix to append to output files. if specified
 * as {} a prefix must be specified on command line */
static s8 g_output_prefix = s8("out/ex");

/* simulation output extent [cm] */
static Rect g_extent = {
	.top   = 1.5,
	.bot   = -1.5,
	.left  = -1.5,
	.right = 1.5,
};

/* incidence location in polar coordinates: r [cm], theta [radians] */
static Vec3Pol g_incidence_location = {
	.r = 1.0,
	.theta = DEG2RAD(30),
};

/* incidence angle between z-axis and the plane normal */
static f64 g_theta_i = DEG2RAD(45); /* [radians] */

/* number of output grid points */
static u32 g_Nx = 64;
static u32 g_Ny = 64;

static u32 g_N_lines            = 8;
static u32 g_N_photons_per_line = 1e6;

/* scattering/absorption coefficients */
static f64 g_mu_a = 0.1;   /* [cm^-1] */
static f64 g_mu_s = 100.0; /* [cm^-1] */

/* scattering anisotropy */
static f64 g_anisotropy = 0.9;

/* refractive indices */
static f64 g_n0 = 1.0;
static f64 g_n  = 1.33;
