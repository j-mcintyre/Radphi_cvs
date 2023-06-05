/*
  $Log: DotProduct3.c,v $
  Revision 1.2  1997/05/29 20:26:40  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: DotProduct3.c,v 1.2 1997/05/29 20:26:40 radphi Exp $";


#include <stdio.h>
#include <math.h>

#include <clib.h>

double DotProduct3(const vector3_t *p1,const vector3_t *p2)
{
    return(p1->x * p2->x + p1->y * p2->y + p1->z * p2->z);
}

/* end file */
