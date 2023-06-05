/* We need algorithm that will sort the hits in CPV,RPD and UPV on the time basis
 * The algorithm for each channel should do the following:
 * find the first hit, consider all 
 * the hits that came within the next 25 ns indistinguishable and save the
 * front edge of the first hit; if the next hit comes after 25 ns
 * it should be considered a separate hit and treated as the first one
 * This subroutine implements the algorithm
*/

#include "digitizer.h"

void timeAnalysis (int *channel,int chdim,int *numbv,float *energy, 
                   int *times,float *timelist,int listdim,float *sorthits,
                   int nhits,int key,int nvdim,int nhdim) {
 int number;
 int i;
 float tstart[100];
 
 if (nhits>0) {
  for (i=0;i<nhits;i++) {
   if (nvdim==1) {
    number=(numbv[0+nvdim*i]-1);
   }
   else { /* In this application it means that nvdim=2; for RPD only */
    number=2*(numbv[0+nvdim*i]-1)+numbv[1+nvdim*i]-1;
   }  
   energy[channel[number+chdim*key]]+=sorthits[4+nhdim*i]; 
   /* accumulate the energy */
   if(times[channel[number+chdim*key]]==0) {
    /*start accumulating hits for this channel*/
    tstart[channel[number+key*chdim]]=sorthits[6+nhdim*i]; 
    /* front edge of the pulse */
    timelist[channel[number+chdim*key]+listdim*times[channel[number+chdim*key]]]=
    sorthits[6+nhdim*i];
    /* for the time-list */
    times[channel[number+chdim*key]]++; /* increase the times of hits saved */
   }
   else { /* this channel has been hit before; check if within the resolution*/
    if ((tstart[channel[number+chdim*key]]+TDC_PULSE_RESOLUTION)<sorthits[6+nhdim*i]) { 
    /* not within the 25 ns */
     timelist[channel[number+chdim*key]+listdim*times[channel[number+chdim*key]]]=
     sorthits[6+nhdim*i]; 
     /* for the time-list, next hit */
     times[channel[number+chdim*key]]++; 
     /* next distinguishable pulse */
     tstart[channel[number+chdim*key]]=sorthits[6+nhdim*i]; 
     /* reset the front edge arrival time */
    }
   }
  }
 }
}

