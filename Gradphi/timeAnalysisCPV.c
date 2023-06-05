/* modified from timeAnalysis.c by Tom Bogue 6/99 */

/* We need algorithm that will sort the hits in CPV on the time basis
 * The algorithm for each channel should do the following:
 * find the first hit, consider all 
 * the hits that came within the next 25 ns indistinguishable and save the
 * front edge of the first hit; if the next hit comes after 25 ns
 * it should be considered a separate hit and treated as the first one
 * This subroutine implements the algorithm
*/

#include <float.h>
#include "digitizer.h"

void timeAnalysisCPV (int *numbv,float *energy, 
                   int *times,float *timelist,int listdim,float *sorthits,
                   int nhits,int key,int nvdim,int nhdim) {
 int number;
 int i;
 float tstart[100];
 float Epulse[100];
 
 for (i=0;i<100;i++) {
  tstart[i]=-999;
 }
 if (nhits>0) {
  for (i=0;i<nhits;i++) {
   number=numbv[i];
   /* for the CPV, the channel number and the volume number are identical */
   energy[number]+=sorthits[4+nhdim*i]; 
   /* accumulate the energy */
   if ((tstart[number]+TDC_PULSE_RESOLUTION)<sorthits[6+nhdim*i]) { 
    /*start accumulating hits for this channel*/
    tstart[number]=sorthits[6+nhdim*i]; 
    Epulse[number]=sorthits[4+nhdim*i]; 
   }
   else {
    Epulse[number]+=sorthits[4+nhdim*i]; 
   }
   if (Epulse[number]>CPV_TDC_THRESHOLD) {
    /* front edge of the pulse */
    timelist[number+listdim*times[number]]=tstart[number];
    /* for the time-list */
    times[number]++; /* increase the times of hits saved */
    /* reset the front edge arrival time */
    Epulse[number]=-FLT_MAX;
   }
  }
 }
}

