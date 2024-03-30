/* global data that will not be modified after startup */
static struct {
	f64 dx, dy;     /* cm */
	u32 Nx, Ny;
	u32 N_photons;

	f64 mu_a, mu_s; /* cm^-1 */
	f64 g, d;

	f64 n, n0;
	f64 theta_i;    /* degrees */

	f64 mu_t;
	f64 xoff, yoff;
} gctx =  {
	.dx = 0.05, .dy = 0.05,
	.Nx = 61, .Ny = 61,
	.N_photons = 10e6,

	.mu_a = 0.1, .mu_s = 100.0,
	.g = 0.9, .d = 1e6,

	.n = 1.0, .n0 = 1.0,
	.theta_i = 45,
};
