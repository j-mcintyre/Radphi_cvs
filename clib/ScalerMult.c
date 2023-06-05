/*
  $Log: ScalerMult.c,v $
  Revision 1.1  1997/05/29 20:25:33  radphi
  Initial revision by lfcrob@dustbunny
  Multiplies a 3 vector by a scaler

  */

static const char rcsid[]="$Id";

#include <stdio.h>
#include <clib.h>

#define D(x) ((double)(x))

void ScalerMult(vector3_t *vout,float s,const vector3_t *vin)
{
  vout->x = vin->x * D(s);
  vout->y = vin->y * D(s);
  vout->z = vin->z * D(s);
}
