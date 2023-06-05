/*
  $Log: flipint.c,v $
  Revision 1.2  1997/05/29 20:26:42  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

  */

#include <clib.h>

static const char rcsid[]="$Id: flipint.c,v 1.2 1997/05/29 20:26:42 radphi Exp $";


void flipint(char *c)

{
	char c1,c2,c3,c4;
	c1 = *c;
	c2 = *(c + 1);
	c3 = *(c + 2);
	c4 = *(c + 3);
	*c = c4;
	*(c + 1) = c3;
	*(c + 2) = c2;
	*(c + 3) = c1;
}
