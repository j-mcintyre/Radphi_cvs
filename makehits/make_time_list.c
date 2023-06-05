/*     make_time_list.c
 *     R.T. Jones    June 9 2000 
 *
 *     Sept. 4, 2003 -RTJ-
 *     Modified to operate on either Monte Carlo or real data.  In the case
 *     of Monte Carlo, correct any bugs in the encoding and save a copy of
 *     the original hits data in case someone needs to look at it.
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <itypes.h> 
#include <iitypes.h>
#include <eventType.h>
#include <ntypes.h> 

#include <disData.h>
#include <pedestals.h>
#include <calibration.h> 
#include <makeHits.h>
#include <triggerType.h>



int make_time_list(itape_header_t *event, time_list_t *time_list, int maxhits)
{
  time_list_t *list;
  int i;

  if (event->trigger == TRIG_MC) {
    list = data_getGroup(event,GROUP_TIME_LIST,0);
    if (list == 0) {
      if (time_list) {
        time_list->nelements = 0;
      }
      return 0;
    }
    else if (event->latch < 200) {	 /* fix any Gradphi encoding bugs */
      time_list_t *mclist;
      mclist = data_getGroup(event,GROUP_TIME_LIST,1);
      if (mclist == 0) {
     /* 
      * Back up the original group to subgroup index 1,
      * then overwrite original group with corrected data.
      * This ensures that the event record contains blessed
      * hits information, even if a copy is returned in time_list.
      */
        int size = sizeof_time_list_t(list->nelements);
        mclist = data_addGroup(event,BUFSIZE,GROUP_TIME_LIST,1,size);
        memcpy(mclist,list,size);
      /*
       * Monte Carlo encoding bug fixes:
       * -times were recorded in seconds, should be ns [rtj]
       */
        for (i=0; i < mclist->nelements; i++) {
          list->element[i].le *= 1e9;
        }
      }
    }
    if (time_list) {
      int size = sizeof_time_list_t(list->nelements);
      memcpy(time_list,list,size);
    }
    return 0;
  }

  if (time_list) {
    time_list->nelements = (p_timelist->nelements < maxhits)?
                            p_timelist->nelements : maxhits;
    for(i=0;i<p_timelist->nelements;i++) {
      time_list->element[i].le = p_timelist->element[i].le;
    }
  }
  else {
    int size = sizeof_time_list_t(p_timelist->nelements);
    time_list_t *tmp = data_addGroup(event,BUFSIZE,GROUP_TIME_LIST,0,size);
    memcpy(tmp,p_timelist,size);
  }

  return(0);
}
