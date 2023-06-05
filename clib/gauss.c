/*
  $Log: gauss.c,v $
  Revision 1.2  1997/05/29 20:26:43  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

  */


#include <math.h>
#include <clib.h>


static const char rcsid[]="$Id: gauss.c,v 1.2 1997/05/29 20:26:43 radphi Exp $";



static float xmean = 0.0,xwidth = 1.0;
static int iset;

void GaussSet(float mean,float width)
{
    xmean = mean;
    xwidth = width;
    iset = 0;
}

float gasdev(void)

{
#define MAXRAN 2147483646

    static float gset;
    float fac,r,v1,v2;
    if (!iset) {
	do {
	    v1 = 2.0 * (float) ( (float) random()/ (float) MAXRAN) - 1.0;
	    v2 = 2.0 * (float) ( (float) random()/ (float) MAXRAN) - 1.0;
	    r = v1 * v1 + v2 * v2;
	}  while (r >= 1.0);
	fac = sqrt(-2.0 * log(r)/r);
	gset = v1 * fac * xwidth + xmean;
	iset = 1;
	return(v2 * fac * xwidth + xmean);
    }
    else {
	iset = 0;
	return(gset);
    }
}
float gauss(float mean,float width)
{

    float fac,r,v1,v2;
	do {
	    v1 = 2.0 * (float) ( (float) random()/ (float) MAXRAN) - 1.0;
	    v2 = 2.0 * (float) ( (float) random()/ (float) MAXRAN) - 1.0;
	    r = v1 * v1 + v2 * v2;
	}  while (r >= 1.0);
	fac = sqrt(-2.0 * log(r)/r);
	return(v2 * fac * width + mean);

}

/* Two-dimensional gaussian random number generator */
void gaussRndm2(double mean1,double mean2,double sigma1,double sigma2,double *rnd1,double *rnd2)
{
  double x1,x2;
  double z1,z2;

  x1 = drand48();
  x2 = drand48();
  z1 = 2.*M_PI*x1;
  z2 = sqrt(-2.*log(x2));
  *rnd1 = sin(z1)*z2 * sigma1 + mean1;
  *rnd2 = cos(z1)*z2 * sigma2 + mean2;
}
