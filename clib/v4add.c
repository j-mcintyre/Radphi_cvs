/*
  $Log: v4add.c,v $
  Revision 1.2  1997/05/29 20:26:48  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v4add.c,v 1.2 1997/05/29 20:26:48 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

vector4_t v4add(const vector4_t *v1,const vector4_t *v2)
{
    vector4_t sum;
    sum.space.x = v1->space.x + v2->space.x;
    sum.space.y = v1->space.y + v2->space.y;
    sum.space.z = v1->space.z + v2->space.z;
    sum.t = v1->t + v2->t;
    return(sum);
}
