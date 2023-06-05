/*
  $Log: p3vector.c,v $
  Revision 1.2  1997/05/29 20:26:45  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: p3vector.c,v 1.2 1997/05/29 20:26:45 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void p3vector(const vector3_t *v,FILE *fp)
{
	fprintf(fp,"x: %lf y: %lf z: %lf\n",v->x,v->y,v->z);
}

