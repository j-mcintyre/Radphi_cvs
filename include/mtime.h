/*
  $Log: mtime.h,v $
  Revision 1.1  1997/03/28 19:19:51  radphi
  Initial revision

  */


#ifndef MTIME_H_INCLUDED
#define MTIME_H_INCLUDED


unsigned long mtime(void);
unsigned long cpumtime(void);
void mwait(unsigned long msec);

#endif
