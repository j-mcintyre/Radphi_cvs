/*
  $Log: flip.c,v $
  Revision 1.2  1997/05/29 20:26:41  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

  */

#include <clib.h>

static const char rcsid[]="$Id: flip.c,v 1.2 1997/05/29 20:26:41 radphi Exp $";

void flip(unsigned char *a, unsigned char *b)
/*
    FUNCTION:
        This routine takes the values of two unsigned characters
        and switches their values.

    INPUT_PARAMETERS:
        a - pointer to unsigned character

        b - pointer to unsigned character
*/
{
	unsigned char temp;
	temp = *b;
	*b = *a;
	*a = temp;
}
