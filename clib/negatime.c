/*
  $Log: negatime.c,v $
  Revision 1.2  1997/05/29 20:26:44  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: negatime.c,v 1.2 1997/05/29 20:26:44 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

vector3_t negative(const vector3_t *v)
{
    vector3_t temp;
    temp.x = -v->x;
    temp.y = -v->y;
    temp.z = -v->z;
    return(temp);
}
