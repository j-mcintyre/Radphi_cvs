/*
  $Log: lorentz.c,v $
  Revision 1.2  1997/05/29 20:26:43  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: lorentz.c,v 1.2 1997/05/29 20:26:43 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

vector4_t lorentz(const vector4_t *beta,const vector4_t *pin)
{ 
    vector4_t ret;
    double d,c,c2;
    d = DotProduct3(&(beta->space),&(pin->space));
    c = beta->t/(beta->t + 1.0);
    c2 = c * d + pin->t;
    ret.space.x = pin->space.x + beta->space.x * beta->t * c2;
    ret.space.y = pin->space.y + beta->space.y * beta->t * c2;
    ret.space.z = pin->space.z + beta->space.z * beta->t * c2;
    ret.t = beta->t * (pin->t + d);
    return(ret);
    
    
}
