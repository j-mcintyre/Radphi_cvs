/*
  $Log: v3sum.c,v $
  Revision 1.2  1997/05/29 20:26:47  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v3sum.c,v 1.2 1997/05/29 20:26:47 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void v3sum(const vector3_t *v1,const vector3_t *v2,vector3_t *sum)
{
	sum->x = v1->x + v2->x;
	sum->y = v1->y + v2->y;
	sum->z = v1->z + v2->z;
}
