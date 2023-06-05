/*
   $Log#
*/

static const char rcsid[] = "$Id: add3.c,v 1.2 1997/05/29 20:26:41 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void add3(vector3_t *v,const vector3_t *plus)
{
	v->x += plus->x;
	v->y += plus->y;
	v->z += plus->z;
}
