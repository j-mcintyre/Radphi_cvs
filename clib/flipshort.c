/*
  $Log: flipshort.c,v $
  Revision 1.2  1997/05/29 20:26:42  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

  */

#include <clib.h>

static const char rcsid[]="$Id: flipshort.c,v 1.2 1997/05/29 20:26:42 radphi Exp $";


void flipshort(char *c)
{
	char c1,c2;
	c1 = *c;
	c2 = *(c + 1);
	*c = c2;
	*(c + 1) = c1;
}
