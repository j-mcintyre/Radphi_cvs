/*
  $Log: v4diff.c,v $
  Revision 1.2  1997/05/29 20:26:48  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v4diff.c,v 1.2 1997/05/29 20:26:48 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void v4diff(const vector4_t *v1,const vector4_t *v2,vector4_t *diff)
{
	v3diff(&v1->space,&v2->space,&diff->space);
	diff->t = v1->t - v2->t;
}
