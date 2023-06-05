/*
  $Log: add4.c,v $
  Revision 1.2  1997/05/29 20:26:41  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: add4.c,v 1.2 1997/05/29 20:26:41 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void add4(vector4_t *v,const vector4_t *plus)
{
	add3(&(v->space),&(plus->space));
	v->t += plus->t;
}
