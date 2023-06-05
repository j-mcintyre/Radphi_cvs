/*
  $Log: mwait.c,v $
  Revision 1.2  1997/05/03 20:38:25  radphi
  Changes made by lfcrob@dustbunny
  Fixed to compile clean under IRIX

 * Revision 1.1  1997/03/28  19:20:42  radphi
 * Initial revision
 *
  */


#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <mtime.h>


void mwait(unsigned long msec)
{
  unsigned int sec;
  unsigned long t1,t2;

  sec = msec/1000;
  
  t1=mtime();
  sleep(sec);
  t2=mtime();
  while((t2-t1) < msec)
    t2=mtime();
  return;
}
