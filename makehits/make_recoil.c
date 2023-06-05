/* make_recoil.c
 *  Craig Steffen  June 20 2001
 *
 * June 2, 2003  R.T. Jones
 * -modified to remove hardwired timing corrections
 *
 * April 21, 2005  R.T. Jones
 * -modified to reflect the following changes in the recoil_t structure:
 *  1) replaced pointer member nrecoilp which was used to indicate which
 *     pixel was supposed to define the t0 for this event with the new
 *     member best_t0 which simply contains the best t0 estimate.
 *  2) slight change to the meaning of "time" member of recoil structure
 *     to be the average time of the three hits making up a pixel, without
 *     any adjustments for average TOF, which were originally put in to
 *     enable the use of "time" as a t0 estimate (see point 1 above).
 *  3) moved the average TOF correction that was formerly applied to the
 *     "time" member to apply it to the best_t0 member.
 *
 * May 17, 2005  R.T. Jones
 *  1) added association of recoil hit to a BGV cluster, requiring
 *     agreement in time and z coordinates.
 *  2) added a new member Etotal to the recoil structure, containing
 *     the sum of BGV and BSD deposited energy (GeV).
 *  3) removed a factor of sin(theta) that was being applied in computing
 *     dEdx from pixel "energy" because it is already inside.
 *
 * August 26, 2005 R.T. Jones
 *  1) added ordering of recoils list in order of increasing time, to
 *     facilitate analysis-time decisions on what TOF cutoff to apply.
 *  2) moved earliest acceptable pixel time (relative to trigger) from
 *     hard-coded value of -50ns to macro EARLIEST_RECOIL_TIME.
 */


#define EARLIEST_RECOIL_TIME -10.0
#define DELTA_T_BGV_PIXEL 3.0
#define DELTA_PHI_BGV_PIXEL 0.200
#define DELTA_Z_BGV_PIXEL 50.0

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <itypes.h> 
#include <iitypes.h>
#include <libBGV.h>
#include <eventType.h>
#include <ntypes.h> 

#include <disData.h>
#include <pedestals.h>
#include <calibration.h> 
#include <makeHits.h>

/* #define DEBUGGING 1 */

#define max(a,b) ((a)>(b))?(a):(b)
#define max3(a,b,c) max(a,max(b,c))
#define nint(a) ((int)(a+0.5))

int make_recoil(itape_header_t *event, recoils_t *p_recoils, int maxrecoils)
{
  recoils_t *recoils;
  bsd_pixels_t *bsd_pixels=NULL;
  bgv_clusters_t *bgv_clusters=NULL;
  int i,j;
  int is,js;
  struct pixel_properties {
    int S;
    int Z;
    int index;
    int cluster_with;
    float mean_time;
    float mean_energy;
    float best_time0;
    float best_energy;
  };
  struct pixel_properties sel_pix[200];
  int N_select_pixels=0;
  float earliest_time;

#ifdef DEBUGGING 
  fprintf(stderr,"make_recoil: run %4d event %6d\n",
          event->runNo,event->eventNo);
#endif

  if (p_recoils) {
    p_recoils->nrecoils = 0;
  }

  bsd_pixels = data_getGroup(event,GROUP_BSD_PIXELS,0);
  if (bsd_pixels == NULL) {
    return(1);
  }

  recoils = malloc(sizeof_recoils_t(maxrecoils));
  recoils->nrecoils = 0;

#ifdef DEBUGGING 
  fprintf(stderr,"make_recoil: looping through %3d BSD pixels:\n",
          bsd_pixels->npixels);
#endif    

  for (i=0; i < bsd_pixels->npixels; i++) {
    float Er,El,Es;
    float Tr,Tl,Ts;
    int S = bsd_pixels->pixel[i].straight - 24;
    int R = bsd_pixels->pixel[i].ring;
#ifdef DEBUGGING 
    fprintf(stderr,"make_recoil: pixel %3d: ",i);
#endif    
    if (R >= BSD_N_RINGS) {
    /* would form pixel off the DS end! */
      fprintf(stderr,"Error in make_recoil: ");
      fprintf(stderr,"found pixel off end of barrel: ");
      return(1);
    }
    if (bsd_pixels->pixel[i].time[0] > 10000.0 ||
        bsd_pixels->pixel[i].time[1] > 10000.0 ||
        bsd_pixels->pixel[i].time[2] > 10000.0) {
#ifdef DEBUGGING 
      fprintf(stderr,"pixel has paddle without tdc hits");
      fprintf(stderr," (happens occasionally)\n");
#endif
      /* this indicates that there is no timing info on at least one paddle
       * in a pixel, so we do not select that pixel.
       */
      continue;
    }
    else if (R < 0 || S < 0 || S >= BSD_N_PIXELS_PER_RING){
      fprintf(stderr,"Error in make_recoil:\n");
      fprintf(stderr,"ring or straight out of range.\n");
      fprintf(stderr,"Z=%2d S=%2d\n",R,S);
      return(1);
    }
    Er = bsd_pixels->pixel[i].energy[0];
    El = bsd_pixels->pixel[i].energy[1];
    Es = bsd_pixels->pixel[i].energy[2];
    Tr = bsd_pixels->pixel[i].time[0];
    Tl = bsd_pixels->pixel[i].time[1];
    Ts = bsd_pixels->pixel[i].time[2];

    if (Ts > EARLIEST_RECOIL_TIME &&         /* if pixel is not too early */
      fabs(Tr-Tl) < 3.0 &&
      fabs(Tl-Ts) < 3.0 &&
      fabs(Ts-Tr) < 3.0) {     /* and if pixel is time correlated */
        sel_pix[N_select_pixels].index = i;
        sel_pix[N_select_pixels].S = S;
        sel_pix[N_select_pixels].Z = R;
        sel_pix[N_select_pixels].mean_time = (Tr+Tl+Ts)/3;
        sel_pix[N_select_pixels].mean_energy = (Er+El+Es)/3;
       /* best_time0 : the best estimate for t0 of this event
        *   The times Tr,Tl,Ts in the pixel structure already have the
        *   light propagation time inside the plastic subtracted away.
        *   The best_time0 is the average of Tr,Tl,Ts with a correction
        *   of (0.5 ns * pixel_ring_index) to subtract average TOF.
        */
        sel_pix[N_select_pixels].best_time0 =
                          (Tr+Tl+Ts)/3. - ((R < 7)? (R-3)*0.5 : 0);
       /* best_energy : the best estimate for dE/dx in the scintillator
        *   The energies Er,El,Es in the pixel structure have already
        *   been corrected for attenuation and path length differences
        *   by dividing by bsd_pixel_gain[].  The best_energy is formed
        *   by discarding the largest of the three dE/dx values and
        *   taking the average of the two remaining values.
        */
        sel_pix[N_select_pixels].best_energy = (Er+El+Es-max3(Er,El,Es))/2;
        N_select_pixels++;
    }
#ifdef DEBUGGING 
    fprintf(stderr,"\n");
#endif
  }

  for (is=0; is < N_select_pixels; is++) {
    sel_pix[is].cluster_with = is;
    for (js=is+1; js < N_select_pixels; js++) {
      float Tdiff = fabs(sel_pix[is].mean_time-sel_pix[js].mean_time);
      if (Tdiff < 3.0) {             /* pixel clustering time window */
        float Zdiff = abs(sel_pix[is].Z-sel_pix[js].Z);
        float Sdiff = abs(sel_pix[is].S-sel_pix[js].S);
        Sdiff = (Sdiff == 23)? 1 : Sdiff;
        if (Zdiff <= 1 && Sdiff <= 1) {
	  sel_pix[is].cluster_with = js;
	}
      }
    }
  }
  for (is=0; is < N_select_pixels; is++) {
    int next = is;
    while (sel_pix[next].cluster_with > next) {
      next = sel_pix[next].cluster_with;
    }
    if (sel_pix[next].cluster_with >= is) {
      int irec = recoils->nrecoils++;
      if (irec >= maxrecoils) {
#ifdef VERBOSE
        fprintf(stderr,"Warning in make_recoil: ");
        fprintf(stderr,"too many recoils for table.\n");
#endif
        irec = --recoils->nrecoils;
      }
      recoils->recoil[irec].npixels = 1;
      recoils->recoil[irec].pixel_index[0] = is;
      sel_pix[next].cluster_with = irec;
    }
    else {
      int irec = sel_pix[next].cluster_with;
      int ipix = recoils->recoil[irec].npixels++;
      if (ipix == MAX_PIXELS_IN_RECOIL) {
#ifdef VERBOSE
        fprintf(stderr,"Warning in make_recoil: ");
        fprintf(stderr,"too many pixels in cluster.\n");
#endif
        recoils->recoil[irec].npixels--;
      }
      else {
        recoils->recoil[irec].pixel_index[ipix] = is;
      }
    }
  }

  earliest_time = 1e9;
  for (i=0; i < recoils->nrecoils; i++) {
    float Z=0;
    float S=0;
    float E=0;
    float T=0;
    float T0=0;
    for (j=0; j < recoils->recoil[i].npixels; j++) {
      is = recoils->recoil[i].pixel_index[j];
      recoils->recoil[i].pixel_index[j] = sel_pix[is].index;
      Z += sel_pix[is].Z;
      S += sel_pix[is].S;
      T += sel_pix[is].mean_time;
      T0 += sel_pix[is].best_time0;
      E += sel_pix[is].best_energy;
    }
    Z /= recoils->recoil[i].npixels;
    S /= recoils->recoil[i].npixels;
    T /= recoils->recoil[i].npixels;
    T0 /= recoils->recoil[i].npixels;
    recoils->recoil[i].theta = atan2(BSD_RADIUS,bsd_ring_midz[nint(Z)]);
    recoils->recoil[i].phi = BSD_PHI_ZERO+(BSD_PHI_INCREMENT*nint(S));
    if (recoils->recoil[i].phi >= TWO_PI)
      recoils->recoil[i].phi -= TWO_PI;
    recoils->recoil[i].time = T;
    recoils->recoil[i].dEdx = E;
    recoils->recoil[i].energy = E*3/sin(recoils->recoil[i].theta);
    if (T0 < earliest_time) {
      earliest_time = T0;
    }
  }
  recoils->best_t0 = (earliest_time < 20.)? earliest_time : 0;


/* now add the BGV cluster information */

  bgv_clusters = data_getGroup(event,GROUP_BGV_CLUSTERS,0);
  if (bgv_clusters) {
    int ic;
    for (ic=0; ic < bgv_clusters->nClusters; ic++) {
      int ir;
      float phi=bgv_clusters->clusters[ic].phi;
      float t=bgv_clusters->clusters[ic].t;
      float z=bgv_clusters->clusters[ic].z;
      for (ir=0; ir < recoils->nrecoils; ir++) {
        float phir=recoils->recoil[ir].phi;
        float dphi=fabs(phi-phir);
        dphi=(dphi<M_PI)?dphi:fabs(2*M_PI-dphi);
        if (dphi < DELTA_PHI_BGV_PIXEL) {
          float tr=recoils->recoil[ir].time;
          tr=(t==0)?t:tr;   // allow BGV clusters without both tdc hits
          if (fabs(t-tr) < DELTA_T_BGV_PIXEL) {
            float zr=BSD_RADIUS/tan(recoils->recoil[ir].theta);
            zr=(t==0)?z:zr;   // allow BGV clusters without both tdc hits
            if (fabs(z-zr) < DELTA_Z_BGV_PIXEL) {
              recoils->recoil[ir].energy += bgv_clusters->clusters[ic].energy;
            }
          }
        }
      }
    }
  }

/* now copy them from the temporary array to output area */
/* and reorder them in order of increasing time */

  if (p_recoils == NULL) {
    int size = sizeof_recoils_t(recoils->nrecoils);
    p_recoils = data_addGroup(event,BUFSIZE,GROUP_RECOIL,0,size);
  }
  p_recoils->nrecoils = 0;
  for (i=0; i < recoils->nrecoils; i++) {
    int isel;
    isel = -1;
    earliest_time = 1e9;
    for (j=0; j < recoils->nrecoils; j++) {
      if (recoils->recoil[j].time < earliest_time) {
        earliest_time = recoils->recoil[j].time;
        isel = j;
      }
    }
    if (isel >= 0) {
      p_recoils->recoil[p_recoils->nrecoils++] = recoils->recoil[isel];
      recoils->recoil[isel].time = 1e9;
    }
    else {
      break;
    }
  }
  p_recoils->best_t0 = recoils->best_t0;
  free(recoils);
  return(0);
}
