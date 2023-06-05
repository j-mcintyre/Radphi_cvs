/*
  $Log: effMass.c,v $
  Revision 1.1  1997/05/29 20:26:16  radphi
  Initial revision by lfcrob@dustbunny
  Gives unsigned effective mass of two particles

  */

static const char rcsid[]="$Id";

#include <stdio.h>
#include <clib.h>
#include <math.h>

float effMass(int npart,const vector4_t **fourvecs)
{
  /* returns 'signed' effective mass */
  vector4_t sum;
  double x;
  float ret;
  sum.space.x = 0.0;
  sum.space.y = 0.0;
  sum.space.z = 0.0;
  sum.t = 0.0;
  while (npart--)
    add4(&sum,fourvecs[npart]);
  x = v4magsq(&sum);
  ret = x > 0 ? sqrt(x) : -sqrt(-x);
  return(ret);
}


