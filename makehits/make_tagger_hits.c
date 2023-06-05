/*     make_tagger_hits.c
 *     R.T. Jones    June 9 2000 
 *
 *     Sept. 4, 2003 -RTJ-
 *     Modified to operate on either Monte Carlo or real data.  In the case
 *     of Monte Carlo, correct any bugs in the encoding and save a copy of
 *     the original hits data in case someone needs to look at it.
 *
 *     Mar. 2, 2004 -RTJ-
 *     Removed any reference to tagger hits in Monte Carlo data.  There is
 *     no tagger in the Gradphi simulation so it makes no sense to work
 *     with hits.  Tagging coincidences in Monte Carlo are faked at the
 *     level of tagger_clusters.
 *
 *     Mar. 24, 2004 -RTJ-
 *     Modified in the case that no ADC or TDC values were found for this
 *     detector to write a hits group with zero entries rather than creating
 *     no group at all, which was the previous (confusing) behavior.
 *
 *     Oct. 25, 2005 -RTJ-
 *     Added an ad-hoc lookup table for time shifts called tagger_shifts.dat
 *     that cannot be placed in the map because the shifts jump around even
 *     within a single run.  They are indexed by run,event.
 */

#define APPLY_TAGGER_SHIFTS 1

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


static int apply_tagger_shifts(itape_header_t *event,
                               tagger_hits_t *p_tag,
                               time_list_t *p_timelist);

int make_tagger_hits(itape_header_t *event, tagger_hits_t *p_tag, int maxhits)
{
  tdc_values_t *tdcs=NULL;
  tagger_hits_t *hits=NULL;
  int hindex[20];
  int i,j;

  if (event->trigger == TRIG_MC) {
    if (p_tag) {
      p_tag->nhits = 0;
    }
    return 0;
  }

  if (p_tag == NULL) {
    hits = malloc(sizeof_tagger_hits_t(maxhits));
  }
  else {
    hits = p_tag;
  }

  if (thisEvent != event->eventNo) {
    if (p_timelist==NULL) {
      p_timelist = malloc(sizeof_time_list_t(10000));
    }
    p_timelist->nelements = 0;
    thisEvent = event->eventNo;
  }

  tdcs = data_getGroup(event,GROUP_TAGGER_TDCS,0);
#ifdef VERBOSE
  if (tdcs == NULL) {
    fprintf(stderr,"Warning in make_tagger_hits: ");
    fprintf(stderr,"No tagger TDC's\n");
  }
#endif

  hits->nhits = 0;
  for (i=0; i<20; hindex[i++]=0) {}

  /* restart random number sequence in same place each time */
  randomNext = 10000*event->runNo + event->eventNo + 92;

  if (tdcs) {
    for (i=0;i<tdcs->ntdc;i++) {
      int end=(tdcs->tdc[i].channel)/19;
      /* end=0 for left tagger counters, end=1 for right counters, 
       * and in the case of the tagger OR (tagger TDC channel 38),
       * end=2
       */
      int channel=(tdcs->tdc[i].channel)%19;
      /* so channel is the index of the bar.  Thus channel 2, end 1 is the 
       * phototube on the opposite end of the same scintillator bar as
       * channel 2, end 0.
       */
      if ((tdcs->tdc[i].le >= TMIN_TAG) &&
	  (tdcs->tdc[i].le <= TMAX_TAG) && 
	  (channel_status_tag[channel] == CHANNEL_STATUS_GOOD) && 
	  (end < 2)) { /* this filters out the tagger OR */
        float value=tdcs->tdc[i].le+RANDOM;
        if (p_timelist->nelements >= MAX_TIME_LIST_LENGTH) {
          fprintf(stderr,"Error in make_tagger_hits: ");
          fprintf(stderr,"Time list max length exceeded!\n");
	  break;
        }
        if (hindex[channel]) {
          int ihit = hindex[channel] - 1;
	  if (hits->hit[ihit].times[end] == 0){
	      hits->hit[ihit].tindex[end] = p_timelist->nelements;
	  }
         /* If we come across a second string of tdc values for a
          * phototube that appeared earlier in the tdc tables, I guess
          * we will orphan the former list and start a new one.
          * Maybe the tdc buffer was not flushed from the last event. -rtj-
          */
	  else if (p_timelist->nelements >
                   hits->hit[ihit].tindex[end] + hits->hit[ihit].times[end]) {
	    hits->hit[ihit].tindex[end] = p_timelist->nelements;
	    hits->hit[ihit].times[end] = 0;
	  }
          p_timelist->element[p_timelist->nelements].le =
	           TDC_1877_GAIN*(value - tag_tdc_offset[tdcs->tdc[i].channel]);
          p_timelist->element[p_timelist->nelements++].le += tag_t0_offset;
	  hits->hit[ihit].times[end]++;
        }
        else if (hits->nhits < maxhits) {
	  hits->hit[hits->nhits].channel = channel;
	  hits->hit[hits->nhits].energy[0] =
	  hits->hit[hits->nhits].energy[1] = tagger_energy_settings[0]
	             * tagger_fraction_by_channel[channel];
	  hits->hit[hits->nhits].times[0] = 0;
	  hits->hit[hits->nhits].times[1] = 0;
	  hits->hit[hits->nhits].times[end] =1;
	  hits->hit[hits->nhits].tindex[end] = p_timelist->nelements;
	  p_timelist->element[p_timelist->nelements].le =
	           TDC_1877_GAIN*(value - tag_tdc_offset[tdcs->tdc[i].channel]); 
	  p_timelist->element[p_timelist->nelements++].le += tag_t0_offset;
	  hindex[channel] = ++(hits->nhits);
        }
        else {
#ifdef VERBOSE
  	  fprintf(stderr,"Warning in make_tagger_hits: ");
	  fprintf(stderr,"buffer overflow at %d hits.\n",maxhits+1);
#endif
	  break;
        }
      }
    }
  }
 
#if APPLY_TAGGER_SHIFTS
  apply_tagger_shifts(event, hits, p_timelist);
#endif

  if (p_tag == NULL) {
    int size = sizeof_tagger_hits_t(hits->nhits);
    tagger_hits_t *tmp = data_addGroup(event,BUFSIZE,GROUP_TAGGER_HITS,0,size);
    memcpy(tmp,hits,size);
    free(hits); 
  }
  return(0);
}

int suppress_dead_tagger_channels(tagger_hits_t *tag)
{
   int i;
   int nok=0;
   int nseen=tag->nhits;
   for (i=0; i<nseen; i++) {
     int channel=tag->hit[i].channel;
     if (channel < 0 || channel >= N_TAG_COUNTERS) {
#ifdef VERBOSE
       fprintf(stderr,"Warning in suppress_dead_tagger_channels:");
       fprintf(stderr,"channel %d out of range.\n",channel);
#endif
     }
     else if (channel_status_tag[channel] == CHANNEL_STATUS_GOOD) {
       if (nok < i) {
         tag->hit[nok].channel = channel;
         tag->hit[nok].energy[0] = tag->hit[i].energy[0];
         tag->hit[nok].energy[1] = tag->hit[i].energy[1];
         tag->hit[nok].times[0] = tag->hit[i].times[0];
         tag->hit[nok].times[1] = tag->hit[i].times[1];
         tag->hit[nok].tindex[0] = tag->hit[i].tindex[0];
         tag->hit[nok].tindex[1] = tag->hit[i].tindex[1];
       }
       ++nok;
     }
   }
   return tag->nhits = nok;
}

static int apply_tagger_shifts(itape_header_t *event,
                               tagger_hits_t *hits,
                               time_list_t *p_timelist)
{
   static FILE *ts = 0;
   static int srunNo, seventNo;
   static float tagger_shift[38];
   int ihit;

   if (ts == 0) {
      char tagger_shifts_file[] = "tagger_shifts.dat";
      char *tagger_shifts_dir = getenv("RADPHI_CONFIG");
      char *tagger_shifts_path;
      if (tagger_shifts_dir == NULL) {
         fprintf(stderr,"apply_tagger_shifts: $RADPHI_CONFIG must be set\n");
         return(1);
      }
      tagger_shifts_path = malloc(strlen(tagger_shifts_dir)+
                                  strlen(tagger_shifts_file)+1);
      sprintf(tagger_shifts_path,"%s/%s",tagger_shifts_dir,tagger_shifts_file);
      if ((ts = fopen(tagger_shifts_path,"r")) == NULL) { 
         fprintf(stderr,"apply_tagger_shifts: file %s not found\n",
                 tagger_shifts_path);
         return(1);
      }
      srunNo = 0;
      seventNo = 0;
   }

   if (event->runNo < srunNo) {
      rewind(ts);
   }
   if (srunNo < event->runNo || seventNo < event->eventNo) {
      char saved_shifts[400];
      char shifts[400];
      char *s;
      int i;
      while (srunNo < event->runNo ||
             (srunNo == event->runNo && seventNo < event->eventNo)) {
         strncpy(saved_shifts,shifts,400);
         if (fscanf(ts,"%d %d %[^\n]\n",&srunNo,&seventNo,shifts) <= 0) {
            fprintf(stderr,"apply_tagger_shifts: run %d not found in db\n",
                    event->runNo);
            return(1);
         }
      }
      if (srunNo != event->runNo) {
         fprintf(stderr,"apply_tagger_shifts: run %d not found in db\n",
                 event->runNo);
         return(1);
      }
      s = (strlen(shifts) > 10)? shifts : saved_shifts;
      for (i=0; i<38; i++) {
         tagger_shift[i] = strtod(s,&s);
      }
   }

   for (ihit=0; ihit<hits->nhits; ihit++) {
      int chanl = hits->hit[ihit].channel;
      int chanr = chanl+19;
      int itl = hits->hit[ihit].tindex[0];
      int itr = hits->hit[ihit].tindex[1];
      int i;
      for (i=0; i<hits->hit[ihit].times[0]; i++) {
         p_timelist->element[itl++].le -= tagger_shift[chanl];
      }
      for (i=0; i<hits->hit[ihit].times[1]; i++) {
         p_timelist->element[itr++].le -= tagger_shift[chanr];
      }
   }
}
