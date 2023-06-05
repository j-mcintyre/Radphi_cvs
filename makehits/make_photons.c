/*     make_photons.c
 *
 *     D.S. Armstrong    March 9 2000
 *     - original version
 *
 *     R.T. Jones	 February 2 2001
 *     - disabled call to correctPhotonEnergy()
 *       because energy, angle correction are now applied at the same time
 *       in makePhotonPatch() -- for more information see TechNote 2001-201.
 *
 *     R.T. Jones	 May 17 2005
 *     - added insertion of barrel neutrals
 *     - require that barrel neutrals lie within some time window of t0
 *       for the event, derived from the best recoil time.
 *     - barrel neutrals must have at least one downstream and one upstream
 *       tdc hit to count as a barrel photon
 */

/* #define OLD_NONLINEAR_CORRECTION    -comment out to enable new corrections */

#define DELTA_T_BGV_NEUTRALS 10.0

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

#ifdef OLD_NONLINEAR_CORRECTION
vector4_t correctPhotonEnergy(vector4_t);
#endif

int setup_makePhotons(itape_header_t *event) {
  int actualRun;
  
#ifdef OLD_NONLINEAR_CORRECTION
  if(ccParameterLoad(PARAMETER_EPSILON,event->runNo,
		       &epsilon,&actualRun)){
    fprintf(stderr,"Error from setup_makePhotons: ");
    fprintf(stderr,"failed to load epsilon parameter from map for run %d\n",
	      event->runNo);
    return(1);
  }
#ifdef DEBUGGING
  fprintf(stderr,"Successfully loaded epsilon from map; real run %d\n",
	    actualRun);
#endif
  if(ccParameterLoad(PARAMETER_BETA,event->runNo,
		       &beta,&actualRun)){
    fprintf(stderr,"Error from setup_makePhotons: ");
    fprintf(stderr,"failed to load beta parameter from map for run %d\n",
	      event->runNo);
    return(1);
  }
#ifdef DEBUGGING
  fprintf(stderr,"Successfully loaded beta from map; real run %d\n",
	    actualRun);
#endif
  }
#endif
  return(0);
}

int make_photons(itape_header_t *event, photons_t *p_phot, int max)
{
  photons_t *photons=NULL;
  vector4_t pGamma[30];
  vector3_t vertex={0.0, 0.0, 1.3};
  int size;
  int i;

  lgd_clusters_t *clusters=data_getGroup(event,GROUP_LGD_CLUSTERS,0);
  bgv_clusters_t *bclusters=data_getGroup(event,GROUP_BGV_CLUSTERS,0);
  recoils_t *recoils=data_getGroup(event,GROUP_RECOIL,0);

  if (p_phot == NULL) {
    photons = malloc(sizeof_photons_t(max));
  }
  else {
    photons = p_phot;
  }
  photons->nPhotons = 0;

  if (clusters == 0) {
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_photons: lgd clusters group not found!");
#endif
  }
  else {   /* form photon 4-vectors for good forward clusters */

    if (clusters->nClusters > max){
#ifdef VERBOSE
      fprintf(stderr,"Warning in make_photons: ");
      fprintf(stderr,"buffer overflow at %d clusters, truncating.\n",max+1);
#endif
      clusters->nClusters = max;
    }

    for (i=0; i<clusters->nClusters; i++) {

  /*  apply corrections to turn clusters into photons   */
      pGamma[i]=makePhotonPatch(&(clusters->cluster[i]),vertex);
 
#ifdef OLD_NONLINEAR_CORRECTION
  /* now apply correction for non-linearity in Pb/glass */
      pGamma[i] = correctPhotonEnergy(pGamma[i]);
#else
  /* energy correction now applied in makePhotonPatch() */
#endif

      photons->photon[i].cluster = i;
      photons->photon[i].energy = pGamma[i].t;
      photons->photon[i].momentum = pGamma[i].space;
      photons->photon[i].position = clusters->cluster[i].space;
      (photons->nPhotons)++;
    }
  }

  if (bclusters == 0) {
#ifdef VERBOSE
    fprintf(stderr,"Error in make_photons: bgv clusters group not found!");
#endif
  }
  else if (recoils == 0) {
#ifdef VERBOSE
    fprintf(stderr,"Error in make_photons: recoil group not found!");
#endif
  }
  else {   /* append 4-vectors for neutral barrel clusters */

    float t0=recoils->best_t0;
    for (i=0; i<bclusters->nClusters; i++) {
      if (bclusters->clusters[i].scintEnergy > 0)
        continue;
      if (photons->nPhotons > max){
#ifdef VERBOSE
        fprintf(stderr,"Warning in make_photons: ");
        fprintf(stderr,"buffer overflow at %d clusters, truncating.\n",max+1);
#endif
      }
      else {
        vector4_t p=makePhotonBGV(&bclusters->clusters[i]);
        float phi=bclusters->clusters[i].phi*(M_PI/180.);
        float t=bclusters->clusters[i].t;
        vector3_t r;
        r.x=43.5*cos(phi);
        r.y=43.5*sin(phi);
        r.z=bclusters->clusters[i].z;
        if ((t != 0) && (fabs(t-t0) < DELTA_T_BGV_NEUTRALS)) {
          photons->photon[photons->nPhotons].cluster = -i-1;
          photons->photon[photons->nPhotons].energy = p.t;
          photons->photon[photons->nPhotons].momentum = p.space;
          photons->photon[photons->nPhotons].position = r;
          (photons->nPhotons)++;
        }
      }
    }
  }

  if (p_phot == NULL) {
    int size = sizeof_photons_t(photons->nPhotons);
    photons_t *tmp = data_addGroup(event,BUFSIZE,GROUP_PHOTONS,0,size);
    memcpy(tmp,photons,size);
    free(photons);
  }

  return(0);
}

#ifdef OLD_NONLINEAR_CORRECTION
vector4_t correctPhotonEnergy(vector4_t photon)
{
    /*   This routine will adjust the reconstructed photon energy to
       account for power-law non-linearity of the LGD response, due
       to attenuation length in the Pb/glass; see Scott Teige's note
       of Dec 15 1999 `Corrections to the Radphi mass scale' 
               D.S. Armstrong
    */

     double factor;
     vector4_t answer;

     if (photon.t > 0) {
       factor = exp((1-epsilon)*log(photon.t) + beta)/photon.t;
       answer.t = factor*photon.t;
       answer.space.x = factor*photon.space.x;
       answer.space.y = factor*photon.space.y;
       answer.space.z = factor*photon.space.z;
       }
    return(answer);
}
#endif
