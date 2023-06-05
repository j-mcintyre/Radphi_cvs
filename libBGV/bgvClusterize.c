#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include <itypes.h>
#include <iitypes.h>
#include <disData.h>
#include <libBGV.h>
#include <makeHits.h>
#include "bgvGeom.h"

/* #define DEBUGGING 1 */

static float Emin=0.06; /* this is the threshold for the cluster seed */
static float Emin2=0.01; /* this is the threshold for outliers */

#define E_UP_WEIGHT 1.
#define E_DOWN_WEIGHT 0.
#define CLUSTER_TIME_WINDOW 10.

void sortzv_(float* a, int* index, int* n,
             int* mode, int* nway,
             int* nsort);

/* finds the time, energy and z of hit hitNum */
void bgvFindtzE(bgv_hit_t *hit,
                time_list_t *time_list,
                float *t, float *z, float *E);

/* finds the phi given a channel number */
float bgvFindPhi(int channelNum);

void bgvSetup(float userEmin,float userEmin2)
{
  Emin=userEmin;
  Emin2=userEmin2;
}

void bgvClusterize(itape_header_t *event,bgv_clusters_t **clusters)
{
  static bgv_cluster_t clust[99]; /* there shouldn't be more than 24 hits */
  int nClusters=0;
  int clusterNum;
  int channelNum;
  int hitNum;
  int hitNum2;
  float energy;
  float t,z,phi;
  int i,i2;

  static float Ehit[99],thit[99],zhit[99];
  static int iord[99];
  int one=1, zero=0;

  bgv_hits_t *hits=data_getGroup(event,GROUP_BGV_HITS,0);
  time_list_t *time_list=data_getGroup(event,GROUP_TIME_LIST,0);
  if (!hits) {
    fprintf(stderr,"Warning in bgvClusterize: ");
    fprintf(stderr,"cannot clusterize BGV without BGV hits group, giving up\n");
    return;
  }
  if (!time_list) {
    fprintf(stderr,"Warning in bgvClusterize: ");
    fprintf(stderr,"cannot clusterize BGV without timing info, giving up.\n");
    return;
  }

  /* sort the hits in order of decreasing energy */
  for (hitNum=0;hitNum<hits->nhits;hitNum++) {
      bgvFindtzE(&hits->hit[hitNum],time_list,
                 &thit[hitNum],&zhit[hitNum],&Ehit[hitNum]);
   }
   sortzv_(Ehit,iord,&hits->nhits,&one,&one,&zero);

  /* loop over all hits */
#ifdef DEBUGGING
  if (hits) {
    fprintf(stderr,"There are %d hits in the BGV hit list:\n",hits->nhits);
  }
  else{
    fprintf(stderr,"The hits structure is empty!\n");
  }
#endif
  for (i=0;i<hits->nhits;i++) {
    hitNum=iord[i]-1;
    if (Ehit[hitNum]>Emin) { /* it's a seed! */
      int highBlock, lowBlock;
      clust[nClusters].t=thit[hitNum];
      clust[nClusters].z=zhit[hitNum];
      clust[nClusters].energy=Ehit[hitNum];
      clust[nClusters].energy_up=hits->hit[hitNum].energy[0];
      clust[nClusters].energy_dn=hits->hit[hitNum].energy[1];
      clust[nClusters].phi=bgvFindPhi(hits->hit[hitNum].channel);
      clust[nClusters].scintEnergy=0;
      /* look for more hits in this cluster */
      highBlock=(hits->hit[hitNum].channel+1)%24;
      lowBlock=(hits->hit[hitNum].channel+23)%24;
      for (i2=i+1;i2<hits->nhits;i2++) {
        hitNum2=iord[i2]-1;
        if ((hits->hit[hitNum2].channel==lowBlock)||
            (hits->hit[hitNum2].channel==highBlock)) {
          energy=Ehit[hitNum2];
          t=thit[hitNum2];
          z=zhit[hitNum2];
          if (energy<Emin2) continue;
          if ((t!=0) && fabs(t-clust[nClusters].t)>2*SIGMA_T_BGV) continue;
          if ((t!=0) && fabs(z-clust[nClusters].z)>2*SIGMA_Z_BGV) continue;
	  /* its a valid part of the cluster so update high/lowBlock */
	  if (hits->hit[hitNum2].channel==highBlock) {
	    highBlock=(highBlock+1)%24;
          }
	  else {
	    lowBlock=(lowBlock+23)%24;
          }
	  phi=bgvFindPhi(hits->hit[hitNum2].channel);
	  /* now add quantities to the cluster */
	  clust[nClusters].phi=(clust[nClusters].energy*clust[nClusters].phi+
			 energy*phi)/(energy+clust[nClusters].energy);
          /* only update z and t if this block has t values! */
	  if (t!=0) {
	      clust[nClusters].z=(clust[nClusters].energy*clust[nClusters].z+
                                energy*z)/(energy+clust[nClusters].energy);
          }
	  if (t!=0) {
	      clust[nClusters].t=(clust[nClusters].energy*clust[nClusters].t+
                                energy*t)/(energy+clust[nClusters].energy);
	  }
	  clust[nClusters].energy+=energy;
	  clust[nClusters].energy_up+=hits->hit[hitNum2].energy[0];
	  clust[nClusters].energy_dn+=hits->hit[hitNum2].energy[1];
	  /* now we must remove this hit from the list of remaining blocks */
	  Ehit[hitNum2]=0;
        }
      } /* hitNum2 for loop */
      nClusters++;
    }
  } /* main for loop */

  if (*clusters==NULL) 
    free(*clusters);
  *clusters=(bgv_clusters_t *)malloc(sizeof_bgv_clusters_t(nClusters));
  (*clusters)->nClusters=nClusters;
  memcpy(&((*clusters)->clusters[0]),clust,sizeof(bgv_cluster_t)*nClusters);
}

void bgvFindtzE(bgv_hit_t *hit,
                time_list_t *time_list,
                float *t,float *z,float *E)
{
  int ntUp=hit->times[0];
  int ntDown=hit->times[1];
  int itUp=hit->tindex[0];
  int itDown=hit->tindex[1];
  int channel=hit->channel;
  float Eup=hit->energy[0];
  float Edown=hit->energy[1];
  int iu,id;
  *z=0;
  *t=0;
  for (id=0;id<ntDown;id++,itDown++) {
    float td=time_list->element[itDown].le;
    if (fabs(td) > 100) continue;
    for (iu=0;iu<ntUp;iu++,++itUp) {
      float tu=time_list->element[itUp].le;
      if (fabs(td-tu) > CLUSTER_TIME_WINDOW) continue;
    /*
     *  The following formula combines time and pulse-height ratio
     *  information to obtain a best estimate for the longtudinal
     *  coordinate of the hit along the length of the calorimeter
     *  module.  The slope of 0.132 was obtained empirically from
     *  a fit to this statistic in pixel-calorimeter coincidences.
     *  The log(Edown/Eup) is not exactly linear in z, but the
     *  magnitude of this term is small (15%) so a nonlinear term
     *  is not justified.
     */
      *z=(tu-td+log(Edown/Eup))/0.132;
      *t=(tu+td)/2;
    }
  }
  Eup*=exp(*z/bgv_atten_length[channel+24]);
  Edown*=exp(-*z/bgv_atten_length[channel]);
  *E=(E_UP_WEIGHT*Eup+E_DOWN_WEIGHT*Edown)/(E_UP_WEIGHT+E_DOWN_WEIGHT);
  if (*E > 1e9) {
    printf("Warning: very high energy in bgvFindtzE (E > 1e9 GeV!)\n");
  }
  *z+=Z_CROSSING_BGV;
}

float bgvFindPhi(int channelNum)
{
  return((float)channelNum*(2*M_PI/24.0));
}
