/* global data that will not be modified after startup */
static struct {
	Rect extent;                 /* extent [cm] */
	Vec3Pol incidence_location;  /* in polar coordinates */
	u32 Nx, Ny;
	u32 N_photons;

	f64 mu_a, mu_s;              /* cm^-1 */
	f64 g, d;

	f64 n, n0;
	f64 theta_i;                 /* radians */

	f64 mu_t;
	f64 xoff, yoff;
	f64 dx, dy;
} gctx =  {
	.extent = (Rect){
		.top   = 1.5,
		.bot   = -1.5,
		.left  = -1.5,
		.right = 1.5
	},
	.incidence_location = (Vec3Pol){ 1.0, DEG2RAD(30) },
	.Nx = 61, .Ny = 61,
	.N_photons = 10e6,

	.mu_a = 0.1, .mu_s = 100.0,
	.g = 0.9, .d = 1e6,

	.n = 1.0, .n0 = 1.0,
	.theta_i = DEG2RAD(45),
};
