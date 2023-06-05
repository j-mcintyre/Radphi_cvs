/*
  $Log: mtime.c,v $
  Revision 1.2  2014/01/06 04:54:08  radphi
   * mtime.c - small changes to get rid of compiler warnings. [rtj]

  Revision 1.1  1997/03/28 19:20:38  radphi
  Initial revision

  */

#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <mtime.h>


static char rcsid[] = "$Id: mtime.c,v 1.2 2014/01/06 04:54:08 radphi Exp $";



unsigned long mtime(void)  /* current time in milliseconds */
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec*1000+tv.tv_usec/1000;
}

unsigned long cpumtime(void) /* user CPU time in milliseconds */            
{
  static long clk_tck = 0;
  static struct tms  tms;

  if (clk_tck==0)
    {
      clk_tck = CLOCKS_PER_SEC;  /* CLOCKS_PER_SEC is a macro the resolves
                                    into a system call, see 'limits.h' */
    }

  times(&tms);

  return (tms.tms_utime * 1000 / clk_tck);
}

/* endfile */
