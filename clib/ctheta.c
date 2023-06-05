/*
  $Log: ctheta.c,v $
  Revision 1.1  1997/05/29 20:25:56  radphi
  Initial revision by lfcrob@dustbunny
  Gives cos(theta) between 2 3-vectors

  */

static const char rcsid[]="$Id: ctheta.c,v 1.1 1997/05/29 20:25:56 radphi Exp $";

#include <stdio.h>
#include <math.h>
#include <clib.h>


float ctheta(const vector3_t*vector1,const vector3_t*vector2)
{
  /* calculate cos of angle between vector1 and vector2
     v1 . v1/(v1 * v2) */
  return(DotProduct3(vector1,vector2)/(v3mag(vector1) * v3mag(vector2)));
}
