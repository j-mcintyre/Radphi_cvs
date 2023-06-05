/*
  $Log: v3sub.c,v $
  Revision 1.2  1997/05/29 20:26:47  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: v3sub.c,v 1.2 1997/05/29 20:26:47 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

vector3_t v3sub(const vector3_t *v1,const vector3_t *v2)
{
  vector3_t sub;
  sub.x = v1->x - v2->x;
  sub.y = v1->y - v2->y;
  sub.z = v1->z - v2->z;
  return(sub);
}

/* end file */
