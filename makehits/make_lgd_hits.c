/*     make_lgd_hits.c
 *     D.S. Armstrong    March 9 2000
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
#include <eventType.h>
#include <ntypes.h> 

#include <disData.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <lgdGeom.h>
#include <pedestals.h>
#include <calibration.h> 
#include <makeHits.h>
#include <detector_params.h>
#include <triggerType.h>

/* #define DEBUGGING 1 */
#define USE_LGD_CC_TUNE 1

    /*   This routine will calculate the hits in the LGD, where
          a hit requires an ADC value sufficiently over threshold,
          and the hit value is in energy units (GeV), i.e. a 
          calibration has been applied.

          The group GROUP_LGD_HITS is then created.

          The routine  setup_makeHits  must be called some time before 
          this routine is invoked, to load geometry and calibrations
          appropriate for the particular run number.
   */

static unsigned char lgd_status_bitmask_internal=LGD_STATUS_BITMASK_DEFAULT;

void set_lgd_status_bitmask(unsigned char bitmask){
  lgd_status_bitmask_internal=bitmask;
}

unsigned char lgd_status_bitmask(void){
  return lgd_status_bitmask_internal;
}

int make_lgd_hits(itape_header_t *event, lgd_hits_t *p_lgd, int maxhits)
{
  
  adc_values_t *lgd=NULL;
  lgd_hits_t *hits=NULL;
  int index;
  unsigned char testmask;

  if (event->trigger == TRIG_MC) {
    hits = data_getGroup(event,GROUP_LGD_HITS,0);
    if (hits == 0) {
      if (p_lgd) {
        p_lgd->nhits = 0;
      }
      return 0;
    }
    else if (event->latch < 0) {	 /* fix any Gradphi encoding bugs */
      lgd_hits_t *mchits;
      mchits = data_getGroup(event,GROUP_LGD_HITS,1);
      if (mchits == 0) {
     /* 
      * Back up the original group to subgroup index 1,
      * then overwrite original group with corrected data.
      * This ensures that the event record contains blessed
      * hits information, even if a copy is returned in p_lgd.
      */
        int size = sizeof_lgd_hits_t(hits->nhits);
        mchits = data_addGroup(event,BUFSIZE,GROUP_LGD_HITS,1,size);
        memcpy(mchits,hits,size);
      /*
       * Monte Carlo encoding bug fixes:
       * -no known bugs as of 4-9-2003 -rtj-
       */
      }
    }
    if (p_lgd) {
      int size = sizeof_lgd_hits_t(hits->nhits);
      memcpy(p_lgd,hits,size);
    }
    return 0;
  }
  if (p_lgd == NULL) {
    hits = malloc(sizeof_lgd_hits_t(maxhits));
  }
  else {
    hits = p_lgd;
  }

  hits->nhits = 0;

  lgd = data_getGroup(event,GROUP_LGD_ADCS,0);
#ifdef VERBOSE
  if (lgd == NULL){
    fprintf(stderr,"Warning in make_lgd_hits: No LGD ADC's\n");
  }
#endif

  /* restart random number sequence in same place each time */
  randomNext = 10000*event->runNo + event->eventNo + 74;

  if (lgd) {
    for (index=0;index<lgd->nadc;index++) {
      if (lgd->adc[index].value > 
         lgd_thresh[lgd->adc[index].channel]) {
#ifdef DEBUGGING 
        fprintf(stderr,"ch %3d val=%#05x\n",
	        lgd->adc[index].channel,lgd->adc[index].value);
        if(lgd->adc[index].channel==219){
	  testmask=UNSIGNED_CHAR_ALL_BITS_SET;
	  testmask |= lgd_status_bitmask_internal;
	  testmask ^= lgd_status_bitmask_internal;
	  fprintf(stderr,
	 	  "ch %3d val=%#05x all_bits=%#04x int_mask=%#04x testmask=%#04x\n",
		  lgd->adc[index].channel,lgd->adc[index].value,
		  UNSIGNED_CHAR_ALL_BITS_SET,
		  lgd_status_bitmask_internal,
		  testmask
		  );
        }
#endif      
        testmask=UNSIGNED_CHAR_ALL_BITS_SET;
        testmask |= lgd_status_bitmask_internal;
        testmask ^= lgd_status_bitmask_internal;

        if(!(channel_status_lgd[lgd->adc[index].channel] & testmask)){
	  if(hits->nhits < maxhits){
	    float value=lgd->adc[index].value+RANDOM;
	    hits->hit[hits->nhits].channel = lgd->adc[index].channel;
	    hits->hit[hits->nhits].energy = 
	      (value - lgd_ped[lgd->adc[index].channel])
#ifdef USE_LGD_CC_TUNE
	      *lgd_cc_tune[lgd->adc[index].channel];
#else
	      *lgd_cc[lgd->adc[index].channel];
#endif
	    (hits->nhits)++;
	  }
	  else{
#ifdef VERBOSE
	    fprintf(stderr,"Warning in make_lgd_hits: ");
	    fprintf(stderr," buffer full at %d hits, truncating.\n",maxhits+1);
#endif
	    break;
	  }
        }
      }
    }
  }
  if (p_lgd == NULL) {
    int size = sizeof_lgd_hits_t(hits->nhits);
    lgd_hits_t *tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_HITS,0,size);
    memcpy(tmp,hits,size);
    free(hits); 
  }

  return(0);
}

int suppress_dead_lgd_channels(lgd_hits_t *lgd)
{
   int i;
   int nok=0;
   int nseen=lgd->nhits;
   for (i=0; i<nseen; i++) {
     int channel=lgd->hit[i].channel;
     if (channel < 0 || channel >= lgd_nChannels) {
#ifdef VERBOSE
       fprintf(stderr,"Warning in suppress_dead_lgd_channels:");
       fprintf(stderr,"channel %d out of range.\n",channel);
#endif
     }
     else if (channel_status_lgd[channel] == CHANNEL_STATUS_GOOD) {
       float Emin=(lgd_thresh[channel] - lgd_ped[channel])
#ifdef USE_LGD_CC_TUNE
	         *lgd_cc_tune[channel];
#else
	         *lgd_cc[channel];
#endif
       if (lgd->hit[i].energy > Emin) {
         if (nok < i) {
           lgd->hit[nok].channel = channel;
           lgd->hit[nok].energy = lgd->hit[i].energy;
         }
         ++nok;
       }
     }
   }
   return lgd->nhits = nok;
}
