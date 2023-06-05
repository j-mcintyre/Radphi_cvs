/*     make_bsd_hits.c
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
 *     must define the macro MC_BSD_TIME_SMEARING below.
 */

#define MC_BSD_TIME_SMEARING 0.60

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

/* #define DEBUGGING 1 */

void rannor_(float *a1,float *a2);

int make_bsd_hits(itape_header_t *event, bsd_hits_t *p_bsd, int maxhits)
{

  adc_values_t *adcs=NULL;
  tdc_values_t *tdcs=NULL;
  bsd_hits_t *hits=NULL;
  int hindex[48];
  int i,j;

  if (event->trigger == TRIG_MC) {
    hits = data_getGroup(event,GROUP_BSD_HITS,0);
    if (hits == 0) {
      if (p_bsd) {
        p_bsd->nhits = 0;
      }
      return 0;
    }
    else if (event->latch < 100) {	 /* BSD channel bug in Gradphi */
      bsd_hits_t *mchits;
      mchits = data_getGroup(event,GROUP_BSD_HITS,1);
      if (mchits == 0) {
     /* 
      * Back up the original group to subgroup index 1,
      * then overwrite original group with corrected data.
      * This ensures that the event record contains blessed
      * hits information, even if a copy is returned in p_bsd.
      */
        int size = sizeof_bsd_hits_t(hits->nhits);
        mchits = data_addGroup(event,BUFSIZE,GROUP_BSD_HITS,1,size);
        memcpy(mchits,hits,size);
      /*
       * Monte Carlo encoding bug fixes:
       * -bsd channel index is off for inner and outer layers -rtj-
       */
	for (i=0; i < hits->nhits; i++) {
          if (hits->hit[i].channel < 12)
            hits->hit[i].channel = (hits->hit[i].channel+1) % 12;
	  else if (hits->hit[i].channel < 24)
            continue;
          else
            hits->hit[i].channel = (hits->hit[i].channel-23) % 24 + 24;
	}
      }
    }
#if defined MC_BSD_TIME_SMEARING
    {	/* introduce channel-independent smearing of simulated data */
      time_list_t *tlist = data_getGroup(event,GROUP_TIME_LIST,0);
      for (i=0; i<hits->nhits; i++) {
        int k=hits->hit[i].tindex;
        for (j=0; j<hits->hit[i].times; j++,k++) {
          float un0, un1;
          rannor_(&un0,&un1);
          tlist->element[k].le += MC_BSD_TIME_SMEARING*un0;
        }
      }
    }
#endif
    if (p_bsd) {
      int size = sizeof_bsd_hits_t(hits->nhits);
      memcpy(p_bsd,hits,size);
    }
    return 0;
  }

  if (p_bsd==NULL) {
    hits = malloc(sizeof_bsd_hits_t(maxhits));
  }
  else {
    hits = p_bsd;
  }

  if (thisEvent != event->eventNo) {
    if (p_timelist==NULL) {
      p_timelist = malloc(sizeof_time_list_t(MAX_TIME_LIST_LENGTH));
    }
    p_timelist->nelements = 0;
    thisEvent = event->eventNo;
  }

  adcs = data_getGroup(event,GROUP_BSD_ADCS,0);
  tdcs = data_getGroup(event,GROUP_BSD_TDCS_LONG,0);
#ifdef VERBOSE
  if ((adcs == NULL) && (tdcs == NULL)) {
    fprintf(stderr,"Warning in make_bsd_hits: ");
    fprintf(stderr,"No bsd ADC's or TDC's\n");
  }
#endif

  hits->nhits = 0;
  for (i=0; i<48; hindex[i++]=0) {}
  
  /* restart random number sequence in same place each time */
  randomNext = 10000*event->runNo + event->eventNo + 4398;

  if (tdcs){
    for (i=0;i<tdcs->ntdc;i++) {
      int channel=tdcs->tdc[i].channel;
      if( (channel < 48) &&      /* makes sure it's an individual channel */
	  (tdcs->tdc[i].le >= TMIN_BSD) &&
	  (tdcs->tdc[i].le <= TMAX_BSD) &&
          (channel_status_bsd[channel]==CHANNEL_STATUS_GOOD)){
        float value=tdcs->tdc[i].le+RANDOM;
        if (p_timelist->nelements >= MAX_TIME_LIST_LENGTH) {
          fprintf(stderr,"Error in make_bsd_hits: ");
          fprintf(stderr,"Time list max length exceeded!\n");
	  break;
        }
        if (hindex[channel]) {
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
#ifdef DEBUGGING
	    fprintf(stderr,"make_bsd_hits: Found duplicate multihit sets ");
            fprintf(stderr,"for channel %d in bsd_tdcs\n", channel);
#endif
          }
          p_timelist->element[p_timelist->nelements].le =
	                  TDC_1877_GAIN*(value - bsd_tdc_offset[channel]);
          p_timelist->element[p_timelist->nelements++].le += bsd_t0_offset;
          hits->hit[ihit].times++;
        }
        else if(hits->nhits < maxhits){
          hits->hit[hits->nhits].channel = channel;
          hits->hit[hits->nhits].energy = 0;
          hits->hit[hits->nhits].times =1;
          hits->hit[hits->nhits].tindex = p_timelist->nelements;
          p_timelist->element[p_timelist->nelements].le =
	                  TDC_1877_GAIN*(value - bsd_tdc_offset[channel]);
          p_timelist->element[p_timelist->nelements++].le += bsd_t0_offset;
          hindex[channel] = ++(hits->nhits);
        }
        else{
#ifdef VERBOSE
          fprintf(stderr,"Warning in make_bsd_hits: ");
          fprintf(stderr,"buffer overflow at %d hits.\n",maxhits+1);
#endif
          break;
        }
      }
    }           /* for (i=0;i<tdcs->ntdc;i++) */
  }             /* if (tdcs) */

  if (adcs) {
    for (i=0;i<adcs->nadc;i++) {
      int channel=adcs->adc[i].channel;
      if ( (adcs->adc[i].value >= AMIN_BSD) && 
           (adcs->adc[i].value >= bsd_thresh[channel]) &&
	   (channel_status_bsd[channel] == CHANNEL_STATUS_GOOD) ) {
        float value=adcs->adc[i].value+RANDOM;
        if (hindex[channel]) {
          int ihit = hindex[channel] - 1;
          hits->hit[ihit].energy = bsd_cc[channel]*(value - bsd_ped[channel]);
        }
        else if(hits->nhits < maxhits){
          hits->hit[hits->nhits].channel = channel;
          hits->hit[hits->nhits].energy = 
                                 bsd_cc[channel]*(value - bsd_ped[channel]);
          hits->hit[hits->nhits].times = 0;
          hits->hit[hits->nhits].tindex = 0;
          hindex[channel] = ++(hits->nhits);
        }
        else{
#ifdef VERBOSE
          fprintf(stderr,"Warning in make_bsd_hits: ");
          fprintf(stderr,"buffer overflow at %d hits.\n",maxhits+1);
#endif
          break;
        }
      }
    }
  }

  if (p_bsd == NULL) {
    int size = sizeof_bsd_hits_t(hits->nhits);
    bsd_hits_t *tmp = data_addGroup(event,BUFSIZE,GROUP_BSD_HITS,0,size);
    memcpy(tmp,hits,size);
    free(hits); 
  }
  return(0);
}

int suppress_dead_bsd_channels(bsd_hits_t *bsd)
{
   int i;
   int nok=0;
   int nseen=bsd->nhits;
   for (i=0; i<nseen; i++) {
     int channel=bsd->hit[i].channel;
     if (channel < 0 || channel >= N_BSD_COUNTERS) {
#ifdef VERBOSE
       fprintf(stderr,"Warning in suppress_dead_bsd_channels:");
       fprintf(stderr,"channel %d out of range.\n",channel);
#endif
     }
     else if (channel_status_bsd[channel] == CHANNEL_STATUS_GOOD) {
       float Emin=(bsd_thresh[channel] - bsd_ped[channel])
	        *bsd_cc[channel];
       if (bsd->hit[i].energy > Emin) {
         if (nok<i) {
           bsd->hit[nok].channel = channel;
           bsd->hit[nok].energy = bsd->hit[i].energy;
           bsd->hit[nok].times = bsd->hit[i].times;
           bsd->hit[nok].tindex = bsd->hit[i].tindex;
         }
         ++nok;
       }
     }
   }
   return bsd->nhits = nok;
}
