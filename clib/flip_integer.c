/*
  $Log: flip_integer.c,v $
  Revision 1.2  1997/05/29 20:26:42  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

  */

#include <clib.h>

static const char rcsid[]="$Id: flip_integer.c,v 1.2 1997/05/29 20:26:42 radphi Exp $";

void flip_integer(unsigned char *integer)
/*
    FUNCTION:
        This routine takes an integer and flips its byte format:
        ---------    ---------
        |1|2|3|4| => |2|1|4|3|
        ---------    ---------

    INPUT_PARAMETERS:
        integer - pointer to unsigned character
            the four byte integer to be flipped
*/
{
	flip(integer,(integer + 3));
	flip((integer + 1) , (integer + 2));
}
