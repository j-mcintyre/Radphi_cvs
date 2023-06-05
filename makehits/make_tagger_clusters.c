/* 2001 June 22 
 *  Craig Steffen
 *
 * Mar. 2, 2004 -Richard Jones
 *  Remove explicit references to the tagger hits from which the clusters
 *  are formed.  For Monte Carlo, generate fake tagger_clusters based on
 *  beam photon information in the MC block.
 */

#define TAGGER_MC_TIME_OFFSET 2.7

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <itypes.h> 
#include <iitypes.h>
#include <eventType.h>
#include <triggerType.h>
#include <ntypes.h> 

#include <disData.h>
#include <pedestals.h>
#include <calibration.h> 
#include <makeHits.h>
#include <mctypes.h>

#define TAGGER_T_WINDOW 80.
#define DELTA_T_CLUSTER 3.0
#define DELTA_T_RL 3.0

/* #define DEBUGGING 1 */

int make_tagger_clusters(itape_header_t *event, tagger_clusters_t *p_clusters,
			 int max_clusters){
  tagger_hits_t *tag_hits=NULL;
  tagger_clusters_t *tag_clusters=NULL;
  time_list_t *time_list=NULL;
  int i,j,k;

#ifdef DEBUGGING
  fprintf(stderr,"make_tagger_clusters: begin \n");
#endif

  if (event->trigger == TRIG_MC) {
    return(generate_tagger_cluster_for_MC_event(event,max_clusters));
  }

  if (p_clusters == NULL) {
    tag_clusters = malloc(sizeof_tagger_clusters_t(max_clusters));
  }
  else{
    tag_clusters = p_clusters;
  }
  if (tag_clusters == NULL){
    fprintf(stderr,"ERROR in make_tagger_clusters: can't allocate memory!\n");
    return(1);
  }

  tag_hits = data_getGroup(event,GROUP_TAGGER_HITS,0);
  if (tag_hits == NULL) {
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_tagger_clusters: ");
    fprintf(stderr,"No tagger hits\n");
#endif
    if (p_clusters == NULL)
      free(tag_clusters);
    return(1);
  }

  if ((time_list = data_getGroup(event,GROUP_TIME_LIST,0)) == NULL){
    fprintf(stderr,"Error in make_tagger_clusters: Failed to get time_list\n");
    return(1);
  }

  tag_clusters->nclusters = 0;
  for (i=0; i < tag_hits->nhits; i++) {
    int ctr = tag_hits->hit[i].channel;
    int iL = tag_hits->hit[i].tindex[0];
    for (j=0; j < tag_hits->hit[i].times[0]; ++j) {
      float tL = time_list->element[iL++].le;
      if (fabs(tL) < TAGGER_T_WINDOW) {		// tagger time outer bounds
        int iR = tag_hits->hit[i].tindex[1];
        for (k=0; k < tag_hits->hit[i].times[1]; ++k) {
          float tR = time_list->element[iR++].le;
#ifdef DEBUGGING
          fprintf(stderr,"make_tagger_clusters: loop j=%2d k=%2d \n",j,k);
#endif
          if (fabs(tL-tR) < DELTA_T_RL) {	// tagger R/L coincidence
	    float tmean = (tL+tR)/2.;
	    int nclu = tag_clusters->nclusters++;
            int ic;
            if (tag_clusters->nclusters > max_clusters) {
#ifdef VERBOSE
	      fprintf(stderr,"Warning in make_tagger_clusters: ");
	      fprintf(stderr,"exceeded max number of tagger clusters!\n");
#endif
	      tag_clusters->nclusters = max_clusters;
	      nclu = max_clusters-1;
	    }
#ifdef DEBUGGING
            fprintf(stderr,"make_tagger_clusters: cluster=%2d\n",
	                                      tag_clusters->nclusters);
#endif
	    tag_clusters->cluster[nclu].channel
		 = tag_hits->hit[i].channel;
	    tag_clusters->cluster[nclu].energy
		 = tagger_fraction_by_channel[ctr]
		   * tagger_energy_settings[0];
	    tag_clusters->cluster[nclu].time
		 = tmean;
            tag_clusters->cluster[nclu].next
	         = 0;
	    for (ic=0; ic < tag_clusters->nclusters-1; ++ic) {
	      float tref=tag_clusters->cluster[ic].time;
	      int ctrref=tag_clusters->cluster[ic].channel;
	      if (ctr == ctrref+1 &&
		  fabs(tmean-tref) < DELTA_T_CLUSTER) {	// cluster coincidence
                tag_clusters->cluster[nclu].next
		   = ic;
		break;
              }
	      else if (ctr == ctrref-1 &&
		  fabs(tmean-tref) < DELTA_T_CLUSTER) {	// cluster coincidence
                tag_clusters->cluster[ic].next
		   = nclu;
		break;
	      }
	    }
	  }
	}
      }
    }
  }

  if (p_clusters == NULL) {
    tagger_clusters_t *tmp;
    int size=sizeof_tagger_clusters_t(tag_clusters->nclusters);
    tmp = data_addGroup(event,BUFSIZE,GROUP_TAGGER_CLUSTERS,0,size);
    memcpy(tmp,tag_clusters,size);
    free(tag_clusters);
  }
#ifdef DEBUGGING
  fprintf(stderr,"make_tagger_clusters: finished, returning 0 \n");
#endif
  return(0);
}

int generate_tagger_cluster_for_MC_event(itape_header_t *event,
					 int max_clusters){
  mc_event_t *mc_event=NULL;
  tagger_clusters_t tagger_clusters;
  float beam_energy;
  float E0,E1,E2;
  int n;

  mc_event = data_getGroup(event,GROUP_MC_EVENT,0);
  if (mc_event == NULL) {
    fprintf(stderr,"HEEEELLLLLP!  generate_tagger_cluster_for_MC_event\n");
    fprintf(stderr,"called on an event where there was no\n");
    fprintf(stderr,"GROUP_MC_EVENT!!!\n");
    return(1);
  }
  if (mc_event->npart != 1) {
#ifdef DEBUGGING
    fprintf(stderr,"HELP!  GROUP_MC_EVENT used to contain just one\n");
    fprintf(stderr,"initial state particle.  This has changed, it now\n");
    fprintf(stderr,"has %d particles.  generate_tagger_cluster_for_MC_event\n");
    fprintf(stderr,"in the file make_tagger_clusters needs to be updated!\n");
#endif
    return(1);
  }
  beam_energy = mc_event->part[0].momentum.t;
#ifdef DEBUGGING
  fprintf(stderr,"make_tagger_clusters: beam energy: %f\n",beam_energy);
#endif

  /* Now, we need to set up the parameters of this fake cluster as 
   * well as we can.
   */

  E0 = tagger_energy_settings[0];
  E1 = E0*tagger_energy_settings[1];
  E2 = E0*tagger_energy_settings[2];
  if (beam_energy>=E1 && beam_energy<=E2) {
    tagger_clusters.nclusters = 1;
    tagger_clusters.cluster[0].time = TAGGER_MC_TIME_OFFSET;
    tagger_clusters.cluster[0].energy = beam_energy;
    for (n=0; n<N_TAGGER_ENERGY_CHANNELS-1; n++) {
      float Ediv = E0*(tagger_fraction_by_channel[n]
                      +tagger_fraction_by_channel[n+1])/2.;
      if (beam_energy > Ediv) {
        break;
      }
    }
    tagger_clusters.cluster[0].channel = n;
    tagger_clusters.cluster[0].next = 0;
#ifdef DEBUGGING
  fprintf(stderr,"generate_tagger_cluster_for_MC_event: \n");
  fprintf(stderr,"generated cluster for single photon.\n");
#endif
  }
  else {
    tagger_clusters.nclusters = 0;
  }
  {
    tagger_clusters_t *tmp;
    int size = sizeof_tagger_clusters_t(tagger_clusters.nclusters);
    tmp = data_addGroup(event,BUFSIZE,GROUP_TAGGER_CLUSTERS,0,size);
    memcpy(tmp,&tagger_clusters,size);
  }

  return(0);
}
