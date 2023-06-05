/*     make_bgv_hits.c
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
 *
 *     Oct. 6, 2005 -RTJ-
 *     Added time smearing for Monte Carlo data.  The Monte Carlo times
 *     were far better than the detector could supply.  To enable it you
 *     must define the macro MC_SMEAR_BGV below.
 */

#define MC_SMEAR_BGV 1

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
#include <libBGV.h>

/* #define DEBUGGING 1 */

/* The following three options are obsolete ways of calculating
 * the time that goes into the time list.  Uncomment any one
 * of the following three lines to revert to old behavior.
 */
//#define NO_TIME_OFFSET 1
//#define USE_MEANDIFF_TIME_OFFSET 1
//#define USE_CPV28_TIME_OFFSET 1

int make_bgv_hits(itape_header_t *event, bgv_hits_t *p_bgv, int maxhits)
{

  adc_values_t *adcs=NULL;
  tdc_values_t *tdcs=NULL;
  bgv_hits_t *hits=NULL;
  int hindex[24];
  int i,j;

  if (event->trigger == TRIG_MC) {
    hits = data_getGroup(event,GROUP_BGV_HITS,0);
    if (hits == 0) {
      if (p_bgv) {
        p_bgv->nhits = 0;
      }
      return 0;
    }
    else if (event->latch < 0) {	 /* fix any Gradphi encoding bugs */
      bgv_hits_t *mchits;
      mchits = data_getGroup(event,GROUP_BGV_HITS,1);
      if (mchits == 0) {
     /* 
      * Back up the original group to subgroup index 1,
      * then overwrite original group with corrected data.
      * This ensures that the event record contains blessed
      * hits information, even if a copy is returned in p_bgv.
      */
        int size = sizeof_bgv_hits_t(hits->nhits);
        mchits = data_addGroup(event,BUFSIZE,GROUP_BGV_HITS,1,size);
        memcpy(mchits,hits,size);
      /*
       * Monte Carlo encoding bug fixes:
       * -no known bugs as of 4-9-2003 -rtj-
       */
      }
    }
#if defined MC_SMEAR_BGV
    {	/* introduce channel-dependent smearing of simulated data */
      time_list_t *tlist = data_getGroup(event,GROUP_TIME_LIST,0);
      MCsmearBGV(hits,tlist);
    }
#endif
    if (p_bgv) {
      int size = sizeof_bgv_hits_t(hits->nhits);
      memcpy(p_bgv,hits,size);
    }
    return 0;
  }

  if (p_bgv == NULL) {
    hits = malloc(sizeof_bgv_hits_t(maxhits));
  }
  else {
    hits = p_bgv;
  }

  if (thisEvent != event->eventNo) {
    if (p_timelist == NULL) {
      p_timelist = malloc(sizeof_time_list_t(10000));
    }
    p_timelist->nelements = 0;
    thisEvent = event->eventNo;
  }

  adcs = data_getGroup(event,GROUP_BGV_ADCS,0);
  tdcs = data_getGroup(event,GROUP_BGV_TDCS,0);
#ifdef VERBOSE
  if ((adcs == NULL) && (tdcs == NULL)){
    fprintf(stderr,"Warning in make_bgv_hits: ");
    fprintf(stderr,"No bgv ADC's or TDC's\n");
  }
#endif

  hits->nhits = 0;
  for (i=0; i<24; hindex[i++]=0) {}

  /* restart random number sequence in same place each time */
  randomNext = 10000*event->runNo + event->eventNo + 498;

  if (tdcs){
    for (i=0;i<tdcs->ntdc;i++) {
      int chanend=tdcs->tdc[i].channel;
      if((tdcs->tdc[i].le >= TMIN_BGV) &&
         (tdcs->tdc[i].le <= TMAX_BGV) &&
         (channel_status_bgv[chanend] == CHANNEL_STATUS_GOOD)) {
        float value=tdcs->tdc[i].le+RANDOM;
        int end=1-chanend/24;
        int channel=chanend%24;
        if (p_timelist->nelements >= MAX_TIME_LIST_LENGTH) {
          fprintf(stderr,"Error in make_bgv_hits: ");
          fprintf(stderr,"Time list max length exceeded!\n");
	  break;
        }
        else if (hindex[channel]) {
          int ihit = hindex[channel] - 1;
          if(hits->hit[ihit].times[end] == 0){
            hits->hit[ihit].tindex[end] = p_timelist->nelements;
          }
         /* If we come across a second string of tdc values for a
          * phototube that appeared earlier in the tdc tables, I guess
          * we will orphan the former list and start a new one.
          * Maybe the tdc buffer was not flushed from the last event. -rtj-
          */
          else if(p_timelist->nelements >
                  hits->hit[ihit].tindex[end] + hits->hit[ihit].times[end]){
            hits->hit[ihit].tindex[end] = p_timelist->nelements;
            hits->hit[ihit].times[end] = 0;
#ifdef DEBUGGING
	    fprintf(stderr,"make_bgv_hits: Found duplicate multihit sets ");
            fprintf(stderr,"for channel %d in bgv_tdcs\n", chanend);
#endif
          }

#if defined NO_TIME_OFFSET
	  p_timelist->element[p_timelist->nelements].le =
                         TDC_1877_GAIN*value;
#elif defined USE_MEANDIFF_TIME_OFFSET
	  if(end == 1){
             p_timelist->element[p_timelist->nelements].le =
	                 TDC_1877_GAIN*(value - bgv_offset_ave[channel]
					      + bgv_mean_diff[channel]/2);
	  }
	  else{
             p_timelist->element[p_timelist->nelements].le =
	                 TDC_1877_GAIN*(value - bgv_offset_ave[channel]
					      - bgv_mean_diff[channel]/2);
	  }
#elif defined USE_CPV28_TIME_OFFSET
	  if(end == 1){
	     p_timelist->element[p_timelist->nelements].le = 
	                 TDC_1877_GAIN*(value - bgv_DS_offset[channel]);
	  }
	  else{
	    p_timelist->element[p_timelist->nelements].le = 
	                 TDC_1877_GAIN*(value - bgv_US_offset[channel]);
	  }
#else
          p_timelist->element[p_timelist->nelements].le =
	                 TDC_1877_GAIN*(value - bgv_tdc_offset[chanend]);
#endif	  
	  p_timelist->element[p_timelist->nelements++].le += bgv_t0_offset;
          hits->hit[ihit].times[end]++;
        }

        else if (hits->nhits < maxhits){ /* first hit of this channel */
          hits->hit[hits->nhits].channel = channel;
          hits->hit[hits->nhits].energy[0] = 0;
          hits->hit[hits->nhits].energy[1] = 0;
          hits->hit[hits->nhits].times[0] =0;
          hits->hit[hits->nhits].times[1] =0;
          hits->hit[hits->nhits].times[end] =1;
          hits->hit[hits->nhits].tindex[end] = p_timelist->nelements;
#if defined NO_TIME_OFFSET
	  p_timelist->element[p_timelist->nelements].le = 
                           TDC_1877_GAIN*value;
#elif defined USE_MEANDIFF_TIME_OFFSET
	  if(end == 1){
             p_timelist->element[p_timelist->nelements].le =
	                 TDC_1877_GAIN*(value - bgv_offset_ave[channel]
					      + bgv_mean_diff[channel]/2);
	  }
	  else{
             p_timelist->element[p_timelist->nelements].le =
	                 TDC_1877_GAIN*(value - bgv_offset_ave[channel]
					      - bgv_mean_diff[channel]/2);
	  }
#elif defined USE_CPV28_TIME_OFFSET
	  if(end == 1){
	     p_timelist->element[p_timelist->nelements].le = 
	                 TDC_1877_GAIN*(value - bgv_DS_offset[channel]);
	  }
	  else{
	    p_timelist->element[p_timelist->nelements].le = 
	                 TDC_1877_GAIN*(value - bgv_US_offset[channel]);
	  }
#else
          p_timelist->element[p_timelist->nelements].le =
	                 TDC_1877_GAIN*(value - bgv_tdc_offset[chanend]);
#endif	  
	  if (value == 0) printf("bingo\n");
	  p_timelist->element[p_timelist->nelements++].le += bgv_t0_offset;
	  hindex[channel] = ++(hits->nhits);
        }

        else{
#ifdef VERBOSE
          fprintf(stderr,"Warning in make_bgv_hits: ");
          fprintf(stderr,"buffer overflow at %d hits.\n",maxhits+1);
#endif
          break;
        }
      }
    }
  } 

  if (adcs) {
    for (i=0;i<adcs->nadc;i++) {
      int chanend=adcs->adc[i].channel;
      if ( (adcs->adc[i].value >= AMIN_BGV) && 
           (adcs->adc[i].value >= bgv_thresh[adcs->adc[i].channel]) &&
           (channel_status_bgv[chanend] == CHANNEL_STATUS_GOOD)) {
        float value=adcs->adc[i].value+RANDOM;
        int end=1-chanend/24;
        int channel=chanend%24;
        if (hindex[channel]) {
          int ihit = hindex[channel] - 1;
          hits->hit[ihit].energy[end] =
                                 bgv_cc[chanend]*(value - bgv_ped[chanend]);
        }
        else if(hits->nhits < maxhits){
          hits->hit[hits->nhits].channel = channel;
          hits->hit[hits->nhits].energy[0] = 0;
          hits->hit[hits->nhits].energy[1] = 0;
          hits->hit[hits->nhits].energy[end] =
                                 bgv_cc[chanend]*(value - bgv_ped[chanend]);
          hits->hit[hits->nhits].times[0] =0;
          hits->hit[hits->nhits].times[1] =0;
          hindex[channel] = ++(hits->nhits);
        }
        else {
#ifdef VERBOSE
          fprintf(stderr,"Warning in make_bgv_hits:");
          fprintf(stderr,"buffer overflow at %d hits.\n",maxhits+1);
#endif
          break;
        }
      }
    }
  }

  if (p_bgv == NULL) {
    int size=sizeof_bgv_hits_t(hits->nhits);
    bgv_hits_t *tmp = data_addGroup(event,BUFSIZE,GROUP_BGV_HITS,0,size);
    memcpy(tmp,hits,size);
    free(hits);
  }
  return(0);
}

int suppress_dead_bgv_channels(bgv_hits_t *bgv)
{
   int i;
   int nok=0;
   int nseen=bgv->nhits;
   for (i=0; i<nseen; i++) {
     int channel=bgv->hit[i].channel;
     if (channel < 0 || channel >= N_BGV_COUNTERS) {
#ifdef VERBOSE
       fprintf(stderr,"Warning in suppress_dead_bgv_channels:");
       fprintf(stderr,"channel %d out of range.\n",channel);
#endif
     }
     else if (channel_status_bgv[channel] == CHANNEL_STATUS_GOOD) {
       float Emin0=(bgv_thresh[channel+24] - bgv_ped[channel+24])
		   *bgv_cc[channel+24];
       float Emin1=(bgv_thresh[channel] - bgv_ped[channel])
                   *bgv_cc[channel];
       if (nok < i) {
         bgv->hit[nok].channel = channel;
         bgv->hit[nok].energy[0] = bgv->hit[i].energy[0];
         bgv->hit[nok].energy[1] = bgv->hit[i].energy[1];
         bgv->hit[nok].times[0] = bgv->hit[i].times[0];
         bgv->hit[nok].times[1] = bgv->hit[i].times[1];
         bgv->hit[nok].tindex[0] = bgv->hit[i].tindex[0];
         bgv->hit[nok].tindex[1] = bgv->hit[i].tindex[1];
       }
       if (bgv->hit[i].energy[0] < Emin0) {
         if (bgv->hit[i].energy[1] > Emin1) {
           bgv->hit[nok].energy[0] = 0;
           ++nok;
	 }
       }
       else if (bgv->hit[i].energy[1] < Emin1) {
         if (bgv->hit[i].energy[0] > Emin0) {
           bgv->hit[nok].energy[1] = 0;
           ++nok;
	 }
       }
       else {
         ++nok;
       }
     }
   }
   return bgv->nhits = nok;
}
