/*
  $Log: v4mag.c,v $
  Revision 1.2  1997/05/29 20:26:48  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v4mag.c,v 1.2 1997/05/29 20:26:48 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

float v4mag(const vector4_t *v)
{
	return( sqrt(v4magsq(v)));
}
