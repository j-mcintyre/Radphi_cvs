/*
 $Log: sum4vectors.c,v $
 Revision 1.2  1997/05/29 20:26:46  radphi
 Changes made by lfcrob@dustbunny
 ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: sum4vectors.c,v 1.2 1997/05/29 20:26:46 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void sum4vectors(int n,const vector4_t *fourvecs[],vector4_t *sum)
{
	sum->space.x = 0.0;
	sum->space.y = 0.0;
	sum->space.z = 0.0;
	sum->t = 0.0;
	while (n--)
	  add4(sum,fourvecs[n]);
}
