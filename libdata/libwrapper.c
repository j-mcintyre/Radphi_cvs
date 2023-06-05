/* 
   $Log: libwrapper.c,v $
   Revision 1.1  1998/06/22 08:06:46  radphi
   Initial revision by jonesrt@jlabs4
   Moved libwrapper.c from tapeWriter to libdata so that it can be
   included in the standard libdata object library - RTJ

 * Revision 1.2  1997/05/08  04:54:50  radphi
 * Changes made by lfcrob@dustbunny
 * Correctly sets event->trigger now
 *
 * Revision 1.1  1997/02/12  15:15:32  lfcrob
 * Initial revision
 *
   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <eventType.h>
#include <triggerType.h>
#include <disData.h>
#include <wrapper.h>

static uint32 defaultRun=0;
static uint32 defaultTime=0;


void wrapperSetDefaults(int type,int value)
{
  switch(type){
  case WRAP_DEF_RUN:
    defaultRun=value;
    break;
  case WRAP_DEF_TIME:
    defaultTime=value;
    break;
  default:
    fprintf(stderr,"wrapperSetDefaults: Unknown default type (%d)\n",type);
    break;
  }
}
int wrapEvent(int *cefBuf,itape_header_t *itape,size_t itapeBufSize)
{
  uint16 eventType;
  uint32 eventNo;
  uint16 trigger;
  int length;
  int *rawGroup;
  size_t size;

  data_newItape(itape);
  getEventInfo(cefBuf,&length,&eventType,&eventNo,&trigger);
  itape->time = defaultTime;
  itape->runNo = defaultRun;
  itape->eventType = eventType;
  switch(itape->eventType){
  case EV_DATA:
    itape->trigger = trigger;
    itape->eventNo = eventNo;
    break;
  case EV_SYNC:
  case EV_PRESTART:
  case EV_GO:
  case EV_PAUSE:
  case EV_END:
  case EV_UNKNOWN:
    itape->trigger=0;
    itape->eventNo=0;
    break;
  default:
    fprintf(stderr,"wrapeEvent: internal error - event type %d\n",eventType);
    return(1);
  }
  /* Now add the CEF event as the raw data */
  size = (length+1)*sizeof(uint32);
  rawGroup = data_addGroup(itape,itapeBufSize,GROUP_RAW,0,size);
  memcpy(rawGroup,cefBuf,size);
  
  return(0);
}

void getEventInfo(int *cefBuf,int *length, uint16 *eventType,uint32  *eventNo,uint16 *trigger)
{
  int *workingBuf = cefBuf;
  bankHeader_t *bh;
  int workingLength;

  *length = *(workingBuf++);
  bh=(void *)(workingBuf++);
  switch(bh->dtype){
  case 0x1:
    /* Some sort of control event probably */
    if(bh->num == 0xcc){
      /* We may be in business...*/
      switch(bh->tag){
      case 16:
	/* SYNC */
	/* Make sure the length is right :) */
	if(*length != 5){
	  *eventType = EV_UNKNOWN;
	  return;
	}
	defaultTime = *(workingBuf++);
	/* All we care about for now */
	*eventType = EV_SYNC;
	return;	
      case 17:
	/* PRESTART */
	/* Make sure the length is right :) */
	if(*length != 4){
	  *eventType = EV_UNKNOWN;
	  return;
	}
	defaultTime = *(workingBuf++);
	defaultRun = *(workingBuf++);
	/* All we care about for now */
	*eventType = EV_PRESTART;
	return;
      case 18:
	/* GO */
	/* Make sure the length is right :) */
	if(*length != 4){
	  *eventType = EV_UNKNOWN;
	  return;
	}
	defaultTime = *(workingBuf++);
	/* All we care about for now */
	*eventType = EV_GO;
	return;	
      case 19:
	/* PAUSE */
	/* Make sure the length is right :) */
	if(*length != 4){
	  *eventType = EV_UNKNOWN;
	  return;
	}
	defaultTime = *(workingBuf++);
	/* All we care about for now */
	*eventType = EV_PAUSE;
	return;	
      case 20:
	/* END */
	/* Make sure the length is right :) */
	if(*length != 4){
	  *eventType = EV_UNKNOWN;
	  return;
	}
	defaultTime = *(workingBuf++);
	/* All we care about for now */
	*eventType = EV_END;
	return;	
      default:
	/* Stuck again */
	*eventType = EV_UNKNOWN;
	return;
      }
    }
    else{
      /* Don't know what it is */
      *eventType = EV_UNKNOWN;
      return;
    }
  case 0x10:
    if(bh->num == 0xcc){
      *trigger = 1 << bh->tag;
      /* Looks like a real event, get the next bank header.
	 First make sure their is another bank, then check
	 that it is an event builder bank. If not, we are stuck.
	 If their is another bank, we know that length must be
	 at least 6 = 1 (first bh) + 5 (Event ID Bank)*/
      if(*length<6){
	/* We're stuck - data w/o EB bank */
	*eventType = EV_UNKNOWN;
	return;
      }
      workingLength = *(workingBuf++); /* Length of 2nd bank */
      if(workingLength !=4){
	/* We're stuck - not an EB bank */
	*eventType = EV_UNKNOWN;
	return;
      }
      bh = (void *)(workingBuf++);
      if((bh->tag == 0xc000) && (bh->dtype == 1) && (bh->num==0)){
	/* Event ID Bank!*/
	*eventNo = *(workingBuf++);
	*eventType = EV_DATA;
	return;
      }
      else{
	/* Stuck again */
	*eventType = EV_UNKNOWN;
	return;
      } 
    }
    else{
      /* Don't know what it is */
      *eventType = EV_UNKNOWN;
      return;
    }
  default:
      /* Don't know what it is */
      *eventType = EV_UNKNOWN;
      return;
  }    
}
