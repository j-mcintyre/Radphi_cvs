/*
  $Log: wrapper.h,v $
  Revision 1.1  1997/02/12 15:12:28  lfcrob
  Initial revision

  */

#ifndef WRAPPER_H_INCLUDED
#define WRAPPER_H_INCLUDED

#include <disData.h>

#define WRAP_DEF_RUN 1
#define WRAP_DEF_TIME 2



typedef struct{
  unsigned short tag;
  unsigned char dtype;
  unsigned char num;
} bankHeader_t;

void wrapperSetDefaults(int type,int value);
int wrapEvent(int *cefBuf,itape_header_t *itape,size_t itapeBufSize);
void getEventInfo(int *cefBuf,int *length,uint16 *eventType,uint32  *eventNo,uint16 *trigger);
#endif
