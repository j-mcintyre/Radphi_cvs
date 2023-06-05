/*
  $Log: make4vector.c,v $
  Revision 1.2  1997/05/29 20:26:43  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: make4vector.c,v 1.2 1997/05/29 20:26:43 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void make4vector(const vector3_t *v3,float mass,vector4_t *v4)
{
	float v2;
	v2 = v3magsq(v3);
	v4->space = *v3;
	v4->t = sqrt(v2 + mass * mass);
}

/* end file */
