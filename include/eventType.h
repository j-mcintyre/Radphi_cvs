/*
  $Log: eventType.h,v $
  Revision 1.2  2000/01/14 18:28:08  radphi
  Changes made by jonesrt@zeus
  added a MCGEN event type for raw event generator files -rtj

 * Revision 1.1  1997/02/12  15:11:31  lfcrob
 * Initial revision
 *
  */


#ifndef EVENTTYPE_H_INCLUDED
#define EVENTTYPE_H_INCLUDED

#define EV_MCGEN 0x0
#define EV_DATA 0x1
#define EV_SYNC 0x2
#define EV_PRESTART 0x4
#define EV_GO 0x8
#define EV_PAUSE 0x10
#define EV_END 0x20
#define EV_UNKNOWN 0x8000
#endif
