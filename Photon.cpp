// Photon.cpp: implementation of the Photon class.
//
//////////////////////////////////////////////////////////////////////

#include "Photon.h"
#include "utils.h"

#include <math.h>
/////////////////////////////////////////


double**  Photon::Rd_xy=NULL;
double	 Photon::Rsp=0.0;
/////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Photon::Photon()
{
	w=1.0;
	nScatter=0;
	dead=false;
	posx=0.0;
	posy=0.0;
	posz=0.0;
////////incident direction
	ux=sin(th_in)*n0/n;
	uy=0.0;
	uz=sqrt(1-ux*ux);
//////////////////////////////////////
	s=0.0;
	db=0.0;
}

Photon::~Photon()
{

}

void Photon::Launch()
{
	double rsp, sinthi, sintht, costhi, costht, raux1, raux2;
	sinthi=sin(th_in);
	costhi=cos(th_in);
	sintht=ux;
	costht=uz;
	raux1=sinthi*costht-costhi*sintht;
	raux1=raux1*raux1;
	raux2=sinthi*costht+costhi*sintht;
	raux2=raux2*raux2;
	rsp=0.5*(raux1/raux2+raux1*(1-raux2)/raux2/(1-raux1));
	Rsp+=w*rsp;
	w=w*(1-rsp);
}

void Photon::CalStep()
{
	if (s < ZERO) {
		double ran = random_uniform();
		s = -log(ran + EPS);
	}
}

void Photon::Caldb()
{
	if (uz<-ZERO)  //travel up to reach upper layer
	{
		db=-posz/uz;
	}
	else if (uz>ZERO) //travel down to reach bottom layer
	{
		db=(d-posz)/uz;
	}
	else
	{
		db=BIGVAL; // A big value
	}
}

bool Photon::HitBoundary()
{
	if(db*mut<=s)
		return true;
	else
		return false;
}

void Photon::Move(double ds) //ds has dimension in cm
{
	posx+=ds*ux;
	posy+=ds*uy;
	posz+=ds*uz;
	s=s-ds*mut;
}

void Photon::Absorb()		// calculate absorption
{
	double dw;
	dw=w*mua/mut;//3.24
	w-=dw;
	if (w<0)
	{
		dead=true;
		w=0.0;
		dw=w+dw; //dw equals the small weight dw=(w-dw)+dw
	}
}

void Photon::Scatter()		// calculate next propagation direction
{
	if (!dead)
	{
		double costheta,fei,ran;
		ran = random_uniform();
		if (g!=0)
		{
			double aa=(1-g*g)/(1-g+2*g*ran);
			costheta=(1+g*g-aa*aa)/(2*g);
			if(costheta<-1) costheta=-1;
			else if (costheta>1) costheta=1;
		}
		else
		{
			costheta=2*ran-1;
		}//(3.28)
		double sintheta,sinfei,cosfei;
		sintheta=sqrt(1-costheta*costheta);
		ran = random_uniform();
		fei=2*PI*ran;//(3.29)
		cosfei=cos(fei);
		if (fei<PI)
		{
			sinfei=sqrt(1-cosfei*cosfei);
		}
		else
		{
			sinfei=-sqrt(1-cosfei*cosfei);
		}
		double uux,uuy,uuz;
		if ((SIGN(uz)*uz)<=(1-ZERO))
		{
			uux=sintheta*(ux*uz*cosfei-uy*sinfei)/sqrt(1-uz*uz)+ux*costheta;
			uuy=sintheta*(uy*uz*cosfei+ux*sinfei)/sqrt(1-uz*uz)+uy*costheta;
			uuz=-sintheta*cosfei*sqrt(1-uz*uz)+uz*costheta;
		}//(3.24)
		else //close to normal propagation
		{
			uux=sintheta*cosfei;
			uuy=sintheta*sinfei;
			uuz=SIGN(uz)*costheta;
		}//(3.30)
		ux=uux;
		uy=uuy;
		uz=uuz;
		nScatter+=1;
	}
}

void Photon::ReflectTransmit() // deal with photon across boundary
{
	double sinalphai,sinalphat,ralphai;
	sinalphai=sqrt(1-(uz)*(uz));
	sinalphat=sinalphai*n/n0; //3.35
	if (sinalphat<ZERO) // close to normal
	{
		ralphai=(n-n0)/(n+n0);
		ralphai=ralphai*ralphai;
	}
	else if (sinalphat>1)// Total reflection
	{
		ralphai=1;
	}
	else
	{
		double A,B;
		A=sinalphai*sqrt(1-sinalphat*sinalphat)-sinalphat*sqrt(1-sinalphai*sinalphai);
		A=A*A;
		B=sinalphai*sqrt(1-sinalphat*sinalphat)+sinalphat*sqrt(1-sinalphai*sinalphai);
		B=B*B;
		ralphai=0.5*(2*A-A*A-A*B)/(B*(1-A));//(3.36)
	}
	double ran = random_uniform();
	if (ran<=ralphai) //rebound to current layer
	{
		uz=-uz;
	}
	else //cross the boundary
	{
		if (uz<-ZERO) //hit upper doundary
		{
			if(posx>-x_offset&&posx<x_offset&&posy>-y_offset&&posy<y_offset)
			{
				Rd_xy[int((posx+x_offset)/dx)][int((posy+y_offset)/dy)]+=w;
			}
			dead=true;
		}
		else //hit bottom boundary
		{
			dead=true;
		}
	}
}

void Photon::Termination()  //Russian roulette survial test
{
	if (!dead)
	{
		double wth=1e-4;
		double m=10;
		if (w<=wth)
		{
			double e = random_uniform();
			if (e>(1/m))
			{
				w=0;
				dead=true;
			}
			else
			{
				w=m*w;
			}
		}
	}
}
