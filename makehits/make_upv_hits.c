/*     make_upv_hits.c
 *     R.T. Jones    June 9 2000 
 *
 *     Sept. 4, 2003 -RTJ-
 *     Modified to operate on either Monte Carlo or real data.  In the case
 *     of Monte Carlo, correct any bugs in the encoding and save a copy of
 *     the original hits data in case someone needs to look at it.
 *
 *     Mar. 24, 2004 -RTJ-
 *     Modified in the case that no ADC or TDC values were found for this
 *     detector to write a hits group with zero entries rather than creating
 *     no group at all, which was the previous (confusing) behavior.
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include <itypes.h> 
#include <iitypes.h>
#include <eventType.h>
#include <ntypes.h> 

#include <disData.h>
#include <pedestals.h>
#include <calibration.h> 
#include <makeHits.h>
#include <detector_params.h>
#include <triggerType.h>


int make_upv_hits(itape_header_t *event, upv_hits_t *p_upv, int maxhits)
{

  adc_values_t *adcs=NULL;
  tdc_values_t *tdcs=NULL;
  upv_hits_t *hits=NULL;
  int hindex[10];
  int i,j;

  if (event->trigger == TRIG_MC) {
    hits = data_getGroup(event,GROUP_UPV_HITS,0);
    if (hits == 0) {
      if (p_upv) {
        p_upv->nhits = 0;
      }
      return 0;
    }
    else if (event->latch < 0) {	 /* fix any Gradphi encoding bugs */
      upv_hits_t *mchits;
      mchits = data_getGroup(event,GROUP_UPV_HITS,1);
      if (mchits == 0) {
     /* 
      * Back up the original group to subgroup index 1,
      * then overwrite original group with corrected data.
      * This ensures that the event record contains blessed
      * hits information, even if a copy is returned in p_upv.
      */
        int size = sizeof_upv_hits_t(hits->nhits);
        mchits = data_addGroup(event,BUFSIZE,GROUP_UPV_HITS,1,size);
        memcpy(mchits,hits,size);
      /*
       * Monte Carlo encoding bug fixes:
       * -no known bugs as of 4-9-2003 -rtj-
       */
      }
    }
    if (p_upv) {
      int size = sizeof_upv_hits_t(hits->nhits);
      memcpy(p_upv,hits,size);
    }
    return 0;
  }
  if(p_upv==NULL){
    hits = malloc(sizeof_upv_hits_t(maxhits));
  }
  else{
    hits = p_upv;
  }

  if (thisEvent != event->eventNo) {
    if (p_timelist == NULL) {
      p_timelist = malloc(sizeof_time_list_t(10000));
    }
    p_timelist->nelements = 0;
    thisEvent = event->eventNo;
  }

  adcs = data_getGroup(event,GROUP_UPV_ADCS,0);
  tdcs = data_getGroup(event,GROUP_UPV_TDCS,0);
#ifdef VERBOSE
  if ((adcs == NULL) && (tdcs == NULL)) {
    fprintf(stderr,"Warning in make_upv_hits: ");
    fprintf(stderr,"No upv ADC's or TDC's\n");
  }
#endif

  hits->nhits = 0;
  for (i=0; i<10; hindex[i++]=0) {}

  /* restart random number sequence in same place each time */
  randomNext = 10000*event->runNo + event->eventNo + 1849;

  if (tdcs){
    for (i=0;i<tdcs->ntdc;i++) {
      int channel = tdcs->tdc[i].channel;
      if((tdcs->tdc[i].le >= TMIN_UPV) &&
         (tdcs->tdc[i].le <= TMAX_UPV) &&
         (channel_status_upv[channel] == CHANNEL_STATUS_GOOD)) {
        float value = tdcs->tdc[i].le+RANDOM;
        if (p_timelist->nelements >= MAX_TIME_LIST_LENGTH) {
          fprintf(stderr,"Error in make_upv_hits: ");
          fprintf(stderr,"Time list max length exceeded!\n");
	  break;
        }
        else if (hindex[channel]) {
          int ihit = hindex[channel] - 1;
          if (hits->hit[ihit].times == 0) {
            hits->hit[ihit].tindex = p_timelist->nelements;
          }
         /* If we come across a second string of tdc values for a
          * phototube that appeared earlier in the tdc tables, I guess
          * we will orphan the former list and start a new one.
          * Maybe the tdc buffer was not flushed from the last event. -rtj-
          */
          else if (p_timelist->nelements >
                   hits->hit[ihit].tindex + hits->hit[ihit].times) {
            hits->hit[ihit].tindex = p_timelist->nelements;
            hits->hit[ihit].times = 0;
          }
          p_timelist->element[p_timelist->nelements].le =
	                     TDC_1877_GAIN*(value - upv_tdc_offset[channel]);
          p_timelist->element[p_timelist->nelements++].le += upv_t0_offset;
          hits->hit[ihit].times++;
        }
        else if (hits->nhits < maxhits){
          hits->hit[hits->nhits].channel = channel;
          hits->hit[hits->nhits].energy = 0;
          hits->hit[hits->nhits].times =1;
          hits->hit[hits->nhits].tindex = p_timelist->nelements;
          p_timelist->element[p_timelist->nelements].le =
	                        TDC_1877_GAIN*(value - upv_tdc_offset[channel]);
          p_timelist->element[p_timelist->nelements++].le += upv_t0_offset;
          hindex[channel] = ++(hits->nhits);
        }
        else {
#ifdef VERBOSE
          fprintf(stderr,"Warning in make_upv_hits: ");
          fprintf(stderr,"buffer overflow at %d hits.\n",maxhits+1);
#endif
          break;
        }
      }
    }
  }

  if (adcs) {
    for (i=0;i<adcs->nadc;i++) {
      int channel = adcs->adc[i].channel;
      if ( (adcs->adc[i].value >= AMIN_UPV) &&
           (adcs->adc[i].value >= upv_thresh[channel]) &&
           (channel_status_upv[channel] == CHANNEL_STATUS_GOOD) ) {
        float value=adcs->adc[i].value+RANDOM;
        if (hindex[channel]) {
          int ihit = hindex[channel] - 1;
          hits->hit[ihit].energy =
                             upv_cc[channel]*(value - upv_ped[channel]);
        }
        else if(hits->nhits < maxhits){
          hits->hit[hits->nhits].channel = channel;
          hits->hit[hits->nhits].energy =
                             upv_cc[channel]*(value - upv_ped[channel]);
          hits->hit[hits->nhits].times = 0;
          hits->hit[hits->nhits].tindex = 0;
          hindex[channel] = ++(hits->nhits);
        }
        else {
#ifdef VERBOSE
          fprintf(stderr,"Warning in make_upv_hits: ");
          fprintf(stderr,"buffer overflow at %d hits.\n",maxhits+1);
#endif
          break;
        }
      }
    }
  }

  if (p_upv == NULL) {
    int size = sizeof_upv_hits_t(hits->nhits);
    upv_hits_t *tmp = data_addGroup(event,BUFSIZE,GROUP_UPV_HITS,0,size);
    memcpy(tmp,hits,size);
    free(hits); 
  }
  return(0);
}

int suppress_dead_upv_channels(upv_hits_t *upv)
{
   int i;
   int nok=0;
   int nseen=upv->nhits;
   for (i=0; i<nseen; i++) {
     int channel=upv->hit[i].channel;
     if (channel < 0 || channel >= N_UPV_COUNTERS) {
#ifdef VERBOSE
       fprintf(stderr,"Warning in suppress_dead_upv_channels:");
       fprintf(stderr,"channel %d out of range.\n",channel);
#endif
     }
     else if (channel_status_upv[channel] == CHANNEL_STATUS_GOOD) {
       float Emin=(upv_thresh[channel] - upv_ped[channel])
	          *upv_cc[channel];
       if (upv->hit[i].energy > Emin) {
         if (nok < i) {
           upv->hit[nok].channel = channel;
           upv->hit[nok].energy = upv->hit[i].energy;
           upv->hit[nok].times = upv->hit[i].times;
           upv->hit[nok].tindex = upv->hit[i].tindex;
         }
         ++nok;
       }
     }
   }
   return upv->nhits = nok;
}
