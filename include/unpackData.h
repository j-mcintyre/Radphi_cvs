/*
  $Log: unpackData.h,v $
  Revision 1.2  2000/07/13 13:25:20  radphi
  Changes made by jonesrt@grendl
  added GROUP_FROZEN for a compacted form of the raw data. -rtj

 * Revision 1.1  1997/05/16  05:39:01  radphi
 * Initial revision by lfcrob@jlabs2
 * Data unpacking user interface
 *
  */

#ifndef UNPACKDATA_H_INCLUDED
#define UNPACKDATA_H_INCLUDED 

#include <disData.h>

int unpackEvent(itape_header_t *event, int bufsize);
void unpackSetLogFile(FILE *fp);
void unpackSetVerbose(int level);
int freezeData(itape_header_t *event, int bufsize);
int thawData(itape_header_t *event, int bufsize);

#endif
