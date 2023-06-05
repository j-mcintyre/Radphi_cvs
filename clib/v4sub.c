/*
   $Log: v4sub.c,v $
   Revision 1.2  1997/05/29 20:26:49  radphi
   Changes made by lfcrob@dustbunny
   ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v4sub.c,v 1.2 1997/05/29 20:26:49 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

vector4_t v4sub(const vector4_t *v1,const vector4_t *v2)
{
    vector4_t sub;
    sub.space.x = v1->space.x - v2->space.x;
    sub.space.y = v1->space.y - v2->space.y;
    sub.space.z = v1->space.z - v2->space.z;
    sub.t = v1->t - v2->t;
    return(sub);
}
