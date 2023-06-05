/*
  $Log: Ctime.c,v $
  Revision 1.2  2014/01/06 04:50:35  radphi
  * Ctime.c [rtj] - small change to remove compiler warning.

  Revision 1.1  1997/05/29 20:25:14  radphi
  Initial revision by lfcrob@dustbunny
  Prints time string w/o trailing \n

  */

static const char rcsid[]="$Id: Ctime.c,v 1.2 2014/01/06 04:50:35 radphi Exp $";

#include <stdio.h>
#include <string.h>
#include <clib.h>
#include <time.h>

char *Ctime(const time_t  *t)
{
  /* ctime without the stupid newline */
  char *g;
  g = (char *)ctime((const time_t *)t);
  *(g + strlen(g) - 1) = 0;
  return(g);
}
