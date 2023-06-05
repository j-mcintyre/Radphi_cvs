/*
  $Log: CrossProduct3.c,v $
  Revision 1.2  1997/05/29 20:26:40  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: CrossProduct3.c,v 1.2 1997/05/29 20:26:40 radphi Exp $";


#include <stdio.h>
#include <math.h>

#include <clib.h>

vector3_t CrossProduct3(const vector3_t *p1,const vector3_t *p2)
{
  vector3_t c;

  c.x = p1->y*p2->z - p1->z*p2->y;
  c.y = -(p1->x*p2->z - p1->z*p2->x);
  c.z = p1->x*p2->y - p1->y*p2->x;

  return c;
}

/* end file */
