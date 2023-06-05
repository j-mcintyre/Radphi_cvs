/*
  $Log: scalerMult4V.c,v $
  Revision 1.2  1997/05/29 20:26:45  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: scalerMult4V.c,v 1.2 1997/05/29 20:26:45 radphi Exp $";

#include <itypes.h>
#include <clib.h>


vector4_t scalerMult4V(vector4_t *vect, float scale)

{
  vector4_t ret={0,0,0,0};

  ret.space.x = vect->space.x*scale;
  ret.space.y = vect->space.y*scale;
  ret.space.z = vect->space.z*scale;
  ret.t = vect->t*scale;
  return(ret);
  
}
