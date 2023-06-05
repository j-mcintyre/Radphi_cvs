/*
  $Log: EffMass.c,v $
  Revision 1.2  1997/05/29 20:26:41  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: EffMass.c,v 1.2 1997/05/29 20:26:41 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

float EffMass(int npart,const vector4_t *fourvecs[])
{
	vector4_t sum;
	sum.space.x = 0.0;
	sum.space.y = 0.0;
	sum.space.z = 0.0;
	sum.t = 0.0;
	while (npart--)
		add4(&sum,fourvecs[npart]);
	return(v4mag(&sum));
}
