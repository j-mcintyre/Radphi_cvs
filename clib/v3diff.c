/*
  $Log: v3diff.c,v $
  Revision 1.2  1997/05/29 20:26:46  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v3diff.c,v 1.2 1997/05/29 20:26:46 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void v3diff(const vector3_t *v1,const vector3_t *v2,vector3_t *diff)
{
	diff->x = v1->x - v2->x;
	diff->y = v1->y - v2->y;
	diff->z = v1->z - v2->z;
}
