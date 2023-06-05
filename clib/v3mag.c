/*
  $Log: v3mag.c,v $
  Revision 1.2  1997/05/29 20:26:47  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id";

#include <stdio.h>
#include <math.h>

#include <clib.h>

#define SQR(x) ((x)*(x))

float v3mag(const vector3_t *p)
{
	return(sqrt(SQR(p->x)+SQR(p->y)+SQR(p->z)));
}
