/*
  $Log: BetaGamma.c,v $
  Revision 1.2  1997/05/29 20:26:40  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: BetaGamma.c,v 1.2 1997/05/29 20:26:40 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

vector4_t BetaGamma(const vector4_t *p)
{
    vector4_t ret;
    double m;
    ret.space.x = p->space.x/p->t;
    ret.space.y = p->space.y/p->t;
    ret.space.z = p->space.z/p->t;
    m = sqrt(v4magsq(p));
    ret.t = p->t/m;
    return(ret);
}
