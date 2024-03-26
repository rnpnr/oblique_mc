/* Providing global variables and functions*/
#include <math.h>
#include <new>
#include <time.h>

//global constants
#ifndef PI
	#define PI 3.14159265
#endif

#ifndef ZERO
	#define ZERO 1.0e-6 //for small deflection angle
#endif

#ifndef n0
	#define n0 1.0 //refractive index of ambient medium
#endif

#ifndef BIGVAL
	#define BIGVAL 1e6//big value
#endif

//global variables
#ifdef GLOBALORIGIN
	#define GLOBAL
#else
	#define GLOBAL extern
#endif
GLOBAL double	dx;
GLOBAL double	dy;
GLOBAL int		Nx;
GLOBAL int		Ny;
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

//generate random data
float	ran3(int *idum);
double	RandomNum(void);


//other auxilliary functions
int min(int d1, int d2);
#define SIGN(x) ((x)>=0?1:-1)


////////allocte dynamic array/////////////////////////
using namespace std;
template <class T>
T* alloc1D(int dim1, T t)
{
	T* p;
	try
	{
		p=new T[dim1];
	}
	catch(bad_alloc)
	{
		cerr<<"Fail to allocate 1D array!"<<endl;
		exit(1);
	}
	for (int tmp1=0;tmp1<dim1;tmp1++)
	{
		p[tmp1]=t;
	}
	return p;
}
template <class T>
void dealloc1D(T* p)
{
	delete [] p;
}
///////////
template <class T>
T** alloc2D(int dim1,int dim2,T t)
{
	T** p;
	try
	{
		p=new T*[dim1];
		for (int tmp1=0;tmp1<dim1;tmp1++)
		{
			p[tmp1]=alloc1D(dim2,t);
		}
	}
	catch (bad_alloc)
	{
		cerr<<"Fail to allocate 2D array!";
	}
	return p;
}
template <class T>
void dealloc2D(int dim1,T** p)
{
	for (int tmp1=0;tmp1<dim1;tmp1++)
	{
		dealloc1D(p[tmp1]);
	}
	delete [] p;
}
/////////
template <class T>
T*** alloc3D(int dim1,int dim2,int dim3,T t)
{
	T*** p;
	try
	{
		p=new T**[dim1];
		for (int tmp1=0;tmp1<dim1;tmp1++)
		{
			p[tmp1]=alloc2D(dim2,dim3,t);
		}
	}
	catch (bad_alloc)
	{
		cerr<<"Fail to allocate 3D array!";
	}
	return p;
}
template <class T>
void dealloc3D(int dim1,int dim2,T*** p)
{
	for (int tmp1=0;tmp1<dim1;tmp1++)
	{
		dealloc2D(dim2,p[tmp1]);
	}
	delete [] p;
}
////////////////////////////////////////////////

