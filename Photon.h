// Photon.h: interface for the Photon class.
// describe the behaviors of photons in Tissue
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PHOTON_H__3C1D58D8_26F1_4E3A_9F95_D5051F1FDE29__INCLUDED_)
#define AFX_PHOTON_H__3C1D58D8_26F1_4E3A_9F95_D5051F1FDE29__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000





class Photon  
{
public:
	double	posx, posy, posz;	//current position of photon
	double	ux, uy, uz;			//current direction of photon
	int		nScatter;			//# of scattering events experienced
	double	w;					//weight
	double	s;					//s is the remaining dimensionless step-size
	double	db;					//the distance from boundary
	bool	dead;				//True if photon is dead, otherwise False
	
	static double**	 Rd_xy;   
	static double	Rsp;			//specular reflectance


	Photon();
	virtual ~Photon();
	void	Launch();	// calculate specular reflectance at the surface
	void	CalStep();  // calculate dimensionless step length
	void	Caldb();	// calculate the travel length to reach the boundary
	bool	HitBoundary(); // True if hitboudary
	void	Move(double ds); //move photon towards boundary
	void	Absorb();		// calculate absorption
	void	Scatter();		// calculate next propagation direction
	void	ReflectTransmit(); // deal with photon across boundary
	void	Termination();  //Russian roulette survial test 
};

#endif // !defined(AFX_PHOTON_H__3C1D58D8_26F1_4E3A_9F95_D5051F1FDE29__INCLUDED_)
