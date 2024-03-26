   /* Hybrid-MC by Li Li, Last modified 03/07
refrence: Wang, "Biomedical optics", chpater3&6 */


/* Programming notes
1. Use Cartesian coordinates.
2. Incident beam in x-z plane.
*/

#include <iostream>
#include <fstream>
#include <string>

#define GLOBALORIGIN
#include "utils.h"
#undef GLOABLORIGIN

#include "Photon.h"


using namespace std;
////////////////////////////////////////////////
int main(void)
{
////////Input Interface////////////////////
	cout<<"*****************************************"<<endl;
	cout<<"Hybrid Monte Carlo for Homogeneous Slab with Oblique incidence"<<endl;
	cout<<"By Li Li, 04/2007"<<endl;
	cout<<"Refrence: Wang, Biomedical optics, chpater3&6"<<endl;
	cout<<"*****************************************"<<endl;
	cout<<endl;

	random_init();

	cout<<"File name(Input:*.mci; Output:*.mco):  ";
	string sFilename, sInfilename, sOutfilename;
	cin>>sFilename;
	sInfilename=sFilename+".mci";
	sOutfilename=sFilename+".mco";
	ifstream Infile;
	Infile.open(sInfilename.c_str());
	if(Infile.fail())
	{
		cerr<<"Error! Can't open input file. Err Code 1"<<endl;
		exit(1);
	}
	//	double P0, R;	//Incident energy/power P0, Beam diameter R
	//	Infile>>P0>>R;
	int nNum_photons;	//# of photons used
	Infile>>nNum_photons;
	Infile>>th_in;
	th_in=th_in/180*PI;
	Infile>>dx>>dy>>Nx>>Ny;
	x_offset=dx*Nx/2;
	y_offset=dy*Ny/2;
	Infile>>n>>mua>>mus>>g>>d;
	mut=mua+mus;
	Infile.close();
/////////////////////////////////////////////

////////for timing///////////////////////////
	time_t start,end;	//simulation beginning and finishing time
	time(&start);
/////////////////////////////////////////////
	int tmp1;

	alloc_mat(&Rd_xy, Nx, Ny);

//////////Propagate photons////////////////////////////////
	for (tmp1=1;tmp1<=nNum_photons;tmp1++)
	{
		Photon photon;
		photon.Launch();
		do
		{
			photon.CalStep();
			photon.Caldb();
			if(photon.HitBoundary())
			{
				photon.Move(photon.db);
				photon.ReflectTransmit();
			}
			else
			{
				photon.Move(photon.s/mut);
				photon.Absorb();
				photon.Scatter();
			}
			photon.Termination();
		} while(!photon.dead);
		if (tmp1%(nNum_photons/10)==0)
		{
			cout<<tmp1<<" photons are done!"<<endl;
		}
	}

	time(&end);
	int CompuTime;
	CompuTime=end-start;
	cout<<"This simulation take "<<CompuTime<<" seconds!"<<endl;
/////////////////////////////////////////////


////////Write data before convolution////////
	ofstream f_out(sOutfilename.c_str());
	f_out<<"/////////////Grid coordinates///////////"<<endl;
	f_out<<"x (cm): ";
	for (tmp1=0;tmp1<Nx;tmp1++)
	{
		f_out<<(tmp1+0.5)*dx-x_offset<<' ';
	}
	f_out<<endl;
	f_out<<"y (cm): ";
	for (tmp1=0;tmp1<Ny;tmp1++)
	{
		f_out<<(tmp1+0.5)*dy-y_offset<<' ';
	}
	f_out<<endl;
	f_out<<"/////////////////////////////////////////"<<endl;


	f_out<<"///////////////Rd_xy/////////////////////"<<endl;
	double scale=nNum_photons*dx*dy;
	double *b = Rd_xy.b;
	for (u32 i = 0; i < Rd_xy.Nx; i++) {
		for (u32 j = 0; j < Rd_xy.Ny; j++)
			f_out << b[j] / scale << ' ';
		f_out << endl;
		b += Rd_xy.Ny;
	}
	f_out<<endl;
	f_out<<"/////////////////////////////////////////"<<endl;
}
