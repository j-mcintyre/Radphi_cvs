/*
  $Log: v4sum.c,v $
  Revision 1.2  1997/05/29 20:26:49  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v4sum.c,v 1.2 1997/05/29 20:26:49 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void v4sum(const vector4_t *v1,const vector4_t *v2,vector4_t *sum)
{
	v3sum(&v1->space,&v2->space,&sum->space);
	sum->t = v1->t + v2->t;
}

