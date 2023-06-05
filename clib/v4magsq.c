/*
  $Log: v4magsq.c,v $
  Revision 1.2  1997/05/29 20:26:48  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v4magsq.c,v 1.2 1997/05/29 20:26:48 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

float v4magsq(const vector4_t *p)
{
    return(p->t * p->t - v3magsq(&p->space));
}
