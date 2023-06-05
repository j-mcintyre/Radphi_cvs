/* 2001 June 22 
   Craig Steffen */

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

/* #define DEBUGGING 1 */

/* #define CPVETO_TIMED_WITH_TAGGER 1 */
#define CPVETO_TIMED_WITH_RECOIL 1

int make_tagger_photons(itape_header_t *event, tagger_photons_t *p_photons,
                        int max_photons)
{
  tagger_photons_t *tagger=NULL;
  tagger_clusters_t *clusters=NULL;
  time_list_t *timelist=NULL;
  cpv_hits_t *cpvhits=NULL;
  recoils_t *recoils=NULL;
  int cpv_recoil_coincidences;
  int cpv_tagger_coincidences;
  int icpv;
  int it;
  
  if (p_photons == 0) {
    tagger = malloc(sizeof_tagger_photons_t(max_photons));
  }
  else {
    tagger = p_photons;
  }

  if ((cpvhits=data_getGroup(event,GROUP_CPV_HITS,0)) == NULL) {
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_tagger_photons: ");
    fprintf(stderr,"failed to get cpv hits group!\n");
#endif
    if (p_photons == 0)
      free(tagger);
    return(1);
  }
  if ((timelist=data_getGroup(event,GROUP_TIME_LIST,0)) == NULL) {
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_tagger_photons: ");
    fprintf(stderr,"failed to get time_list for make_tagger_photons\n");
#endif
    if (p_photons == 0)
      free(tagger);
    return(1);
  }
  if ((recoils=data_getGroup(event,GROUP_RECOIL,0)) == NULL) {
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_tagger_photons: ");
    fprintf(stderr,"failed to get recoil group!\n");
#endif
    if (p_photons == 0)
      free(tagger);
    return(1);
  }
  if ((clusters=data_getGroup(event,GROUP_TAGGER_CLUSTERS,0)) == NULL) {
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_tagger_photons: ");
    fprintf(stderr,"failed to get tagger clusters group!\n");
#endif
    if (p_photons == 0)
      free(tagger);
    return(1);
  }
  if (clusters->nclusters > max_photons) {
    fprintf(stderr,"Error in make_tagger_photons: ");
    fprintf(stderr,"too many tagger clusters (%d)\n",
	    clusters->nclusters);
    if (p_photons == 0)
      free(tagger);
    return(1);
  }

  cpv_recoil_coincidences = 0;
  for (icpv=0; icpv < cpvhits->nhits; icpv++) {
    int ii, itl=cpvhits->hit[icpv].tindex;
    for (ii=0; ii < cpvhits->hit[icpv].times; ii++) {
      float dt = timelist->element[itl++].le - recoils->best_t0;
      if (dt > -3. && dt < 3.) {            // cpv-recoil coincidence 
	++cpv_recoil_coincidences;
      }
    }
  }

  tagger->nphotons = 0;
  for (it=0; it < clusters->nclusters; it++) {
    if (clusters->cluster[it].next == 0) {
      int phot = tagger->nphotons++;
      int chan = clusters->cluster[it].channel;
      float tagt = clusters->cluster[it].time;
      float tagE = clusters->cluster[it].energy;
      int coincidences = 0;
      int accidentals = 0;
      float dt = recoils->best_t0-tagt;
      if (dt > -1. && dt < 5.) {              // recoil-tagger coincidence 
        ++coincidences;
      }
      else if (dt > 5. && dt < 11.) {         // recoil-tagger accidentals
        ++accidentals;
      }
      cpv_tagger_coincidences = 0;
      for (icpv=0; icpv < cpvhits->nhits; icpv++) {
        int ii, itl=cpvhits->hit[icpv].tindex;
        for (ii=0; ii < cpvhits->hit[icpv].times; ii++) {
          float dttag = timelist->element[itl++].le - tagt;
          if (dttag > -1. && dttag < 5.) {     // cpv-tagger coincidence 
	    ++cpv_tagger_coincidences;
          }
        }
      }
      tagger->photon[phot].channel = chan;
      tagger->photon[phot].time = tagt;
      tagger->photon[phot].energy = tagE;
      tagger->photon[phot].tagging_weight = (coincidences-accidentals)*(
#if CPVETO_TIMED_WITH_RECOIL
            (cpv_recoil_coincidences == 0)?
#elif CPVETO_TIMED_WITH_TAGGER
            (cpv_tagger_coincidences == 0)? 
#else
#error "You need to chose CPVETO_TIMED_WITH_RECOIL or CPV_TIMED_WITH_TAGGER, or else define some third option."
#endif
            +1 : -tagger_CPleakage_by_channel[chan]);
    }
  }

  if (p_photons == NULL) {
    int size=sizeof_tagger_photons_t(tagger->nphotons);
    tagger_photons_t *tmp;
    tmp = data_addGroup(event,BUFSIZE,GROUP_TAGGER_PHOTONS,0,size);
  
#ifdef DEBUGGING
  fprintf(stderr,"make_tagger_photons: about to add GROUP_TAGGER_PHOTONS\n");
  fprintf(stderr,"with %d photons and class=%d.\n",
	  tagger_photons->nphotons,
	  tagger_photons->class);
#endif

    memcpy(tmp,tagger,size);
    free(tagger);
  }
  return(0);
}
