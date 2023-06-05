/*
  $Log: flip_short.c,v $
  Revision 1.2  1997/05/29 20:26:42  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

  */

#include <clib.h>

static const char rcsid[]="$Id: flip_short.c,v 1.2 1997/05/29 20:26:42 radphi Exp $";

void flip_short(unsigned char *short_int)
/*
    FUNCTION:
        This routine takes a short and flips its byte format:
        -----    -----
        |1|2| => |2|1|
        -----    -----

    INPUT_PARAMETERS:
        short_int - pointer to unsigned character
            the four byte short to be flipped
*/
{
	flip(short_int,(short_int +1));
}
