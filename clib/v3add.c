/*
  $Log: v3add.c,v $
  Revision 1.2  1997/05/29 20:26:46  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v3add.c,v 1.2 1997/05/29 20:26:46 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

vector3_t v3add(const vector3_t *v1,const vector3_t *v2)
{
    vector3_t sum;
    sum.x = v1->x + v2->x;
    sum.y = v1->y + v2->y;
    sum.z = v1->z + v2->z;
    return(sum);
}
