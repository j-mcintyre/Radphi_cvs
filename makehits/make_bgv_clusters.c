/*     make_bgv_clusters.c            */
/*     R.T. Jones      April 22, 2005 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <math.h>

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
#include <libBGV.h>

#define DELTA_T_BGV_BSD 5.0
#define DELTA_PHI_BGV_BSD 0.15

#define SQR(x) ((x)*(x))

int make_bgv_clusters(itape_header_t *event, bgv_clusters_t *p_clust, int max)
{
  bgv_clusters_t *clusters=NULL;
  itape_header_t *tmp;

  bgv_hits_t *bgvHits = data_getGroup(event,GROUP_BGV_HITS,0);
  bsd_hits_t *bsdHits = data_getGroup(event,GROUP_BSD_HITS,0);
  time_list_t *timelist = data_getGroup(event,GROUP_TIME_LIST,0);

/*  associate hits into clusters  */
  if (bgvHits) {
    bgvClusterize(event,&clusters);
  }
  else {
    static int first_time=1;
    if (first_time) {
#ifdef VERBOSE
      fprintf(stderr,"Warning from make_bgv_clusters: ");
      fprintf(stderr,"cannot clusterize without bgv hits.\n");
#endif
      first_time=0;
    }
    clusters=malloc(sizeof_bgv_clusters_t(0));
    clusters->nClusters=0;
  }

/*  associate bgv clusters to bsd hits  */
  if (bsdHits == 0) {
    static int first_time=1;
    if (first_time) {
#ifdef VERBOSE
      fprintf(stderr,"Warning from make_bgv_clusters: ");
      fprintf(stderr,"cannot clusterize without bsd hits.\n");
#endif
      first_time=0;
    }
  }
  else {
    int ic;
    for (ic=0; ic < clusters->nClusters; ic++) {
      int is;
      float phi=clusters->clusters[ic].phi;
      float t=clusters->clusters[ic].t;
      float z=clusters->clusters[ic].z;
      clusters->clusters[ic].scintEnergy=0;
      for (is=0; is < bsdHits->nhits; is++) {
        int sc=bsdHits->hit[is].channel-24;
        if (sc >= 0) {
          float phis=((sc+1)%24+0.5)*(2*M_PI/24);
          float dphi=fabs(phi-phis);
          dphi=(dphi<M_PI)?dphi:fabs(2*M_PI-dphi);
          if (dphi < DELTA_PHI_BGV_BSD) {
            int pt=bsdHits->hit[is].tindex;
            int it;
            for (it=0; it < bsdHits->hit[is].times; it++) {
              if (fabs(timelist->element[pt++].le-t) < DELTA_T_BGV_BSD) {
                clusters->clusters[ic].scintEnergy+=bsdHits->hit[is].energy;
              }
            }
          }
        }
      }
    }
  }

  if (p_clust == NULL) {
    int size = sizeof_bgv_clusters_t(clusters->nClusters);
    bgv_clusters_t *tmp = data_addGroup(event,BUFSIZE,GROUP_BGV_CLUSTERS,0,size);
    memcpy(tmp,clusters,size);
  }
  else if (clusters->nClusters <= max) {
    memcpy(p_clust,clusters,sizeof_bgv_clusters_t(clusters->nClusters));
  }
  else {
    memcpy(p_clust,clusters,sizeof_bgv_clusters_t(max));
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_bgv_clusters: ");
    fprintf(stderr,"buffer overflow with %d clusters, truncating.\n",
                      clusters->nClusters);
#endif
  }

  free(clusters); 
  return(0);
}
