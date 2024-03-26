// Photon.h: interface for the Photon class.
// describe the behaviors of photons in Tissue
//////////////////////////////////////////////////////////////////////

#ifndef _PHOTON_H
#define _PHOTON_H

class Photon
{
public:
	double posx, posy, posz; // current position of photon
	double ux, uy, uz;       // current direction of photon
	int nScatter;            // # of scattering events experienced
	double w;                // weight
	double s;                // s is the remaining dimensionless step-size
	double db;               // the distance from boundary
	bool dead;               // True if photon is dead, otherwise False

	static double Rsp;       // specular reflectance

	Photon();
	virtual ~Photon();
	void Launch();          // calculate specular reflectance at the surface
	void CalStep();         // calculate dimensionless step length
	void Caldb();           // calculate the travel length to reach the boundary
	bool HitBoundary();     // True if hitboudary
	void Move(double ds);   // move photon towards boundary
	void Absorb();          // calculate absorption
	void Scatter();         // calculate next propagation direction
	void ReflectTransmit(); // deal with photon across boundary
	void Termination();     // Russian roulette survial test
};

#endif // _PHOTON_H
