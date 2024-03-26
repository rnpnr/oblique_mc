#include "utils.h"

/***********************************************************
* A random number generator from Numerical Recipes in C. ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9
float ran3(int &idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	if (idum < 0 || iff == 0)
	{
		iff=1;
		mj=MSEED-(idum < 0 ? -idum : idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++)
		{
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++)
			{
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		idum=1;
	}
	if (++inext == 56)
		inext=1;
	if (++inextp == 56)
		inextp=1;
		mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
		ma[inext]=mj;
		return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
 *	Generate a random number between 0 and 1.  Take a
 *	number as seed the first time entering the function.
 *	The seed is limited to 1<<15.
 *	Wang et al found that when idum is too large, ran3 may return
 *	numbers beyond 0 and 1.
 ****/
//#define STANDARDTEST 0  /* testing program using fixed rnd seed. */
double RandomNum(void)
{
  static bool first_time=1;
  static int idum;	/* seed for ran3. */
  if(first_time)
  {
	#if STANDARDTEST /* Use fixed seed to test the program. */
		idum = - 1;
	#else
		idum = -(int)time(NULL)%(1<<15); /* use 16-bit integer as the seed. */
	#endif
    ran3(idum);
    first_time = 0;
    idum = 1;
  }
  return( (double)ran3(idum) );
}

/////////////////////////////////////////////
int min(int d1,int d2)
{
	if (d1<d2)
	{
		return d1;
	}
	else
	{
		return d2;
	}
}
