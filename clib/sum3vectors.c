/*
  $Log: sum3vectors.c,v $
  Revision 1.2  1997/05/29 20:26:46  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/



#include <stdio.h>
#include <math.h>

#include <clib.h>

static const char rcsid[]="$Id: sum3vectors.c,v 1.2 1997/05/29 20:26:46 radphi Exp $";

void sum3vectors(int n,vector3_t *threevecs[],vector3_t *sum)
{
  
  sum->x = 0.0; sum->y = 0.0; sum->z = 0.0;
     
	while (n--)
	  add3(sum,threevecs[n]);
}
