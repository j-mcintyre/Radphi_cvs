/* We need a special algorithm to sort the BSD hits, since they also need to 
 * be broken up into individual scintillators.
 * The only change from timeAnalysis is the routine for finding number.
 * The algorithm for each channel should do the following:
 *   -find the first hit
 *   - consider all the hits that came within the next 25 ns indistinguishable
 *     and save the front edge of the first hit
 *   - if the next hit comes after 25 ns it should be considered a separate
 *     hit and treated as the first one
 * This subroutine implements the algorithm.
*/

#include "BSDgeom.h"
#include "digitizer.h"
#include <float.h>
#include <math.h>

void timeAnalysisBSD (int *channel,int chdim,int *numbv,float *energy, 
                   int *times,float *timelist,int listdim,float *sorthits,
                   int nhits,int key,int nvdim,int nhdim) {
 int number;
 int i;
 float tstart[100];
 float Epulse[100];
/* the ring number (2 for BSR2, 3 for BSR3 (since its the outer ring))*/
 int RingNum;
/* the scintillator number (starts at 1 for each ring) */
 int ScintNum;
/* the starting angle of the scint */
 float phibase;
/* the current angle of the scint */
 float phi;
/* the coodinates of the hit */
 float x,y,z;
/* the number of scints in this ring */
 int NumScints;
/* the twist of this ring */
 float Twist;
/* the initial angle of this ring */
 float phi0;

 for (i=0;i<100;i++) {
   tstart[i]=-999;
 }
 if (nhits>0) {
  for (i=0;i<nhits;i++) {
/* number = (ring # - 1)*24 + scint. # - 1
 * i.e. for the 6th scint of BSR2, number = (2-1)*24+6-1=29
 */
   RingNum=key+1;
/* key is set in BSDsave to be the ring number - 1*/
   x=sorthits[0+i*nhdim];
   y=sorthits[1+i*nhdim];
   z=sorthits[2+i*nhdim];
   if (fabs(x)<1e-6)
     phi = (y<0) ? -3.1415926535897932/2.0 : +3.1415926535897932/2.0;
   else
     phi=atan2(y,x);
   phi*=180.0/3.1415926535897932;
   if (y<0.0) phi+=360.0;
/* phi = phi0 + Twist*(z-z0)/zmax  
 * this implies that
 * phi0 = phi - Twist*(z-z0)/zmax
 */
   switch(RingNum)
     {
     case 1 : 
       NumScints=NumScintsRing1;
       Twist=BSDtwistRing1;
       phi0=BSDphi0Ring1;
       break;
     case 2 :
       NumScints=NumScintsRing2;
       Twist=BSDtwistRing2;
       phi0=BSDphi0Ring2;
       break;
     case 3 :
       NumScints=NumScintsRing3;
       Twist=BSDtwistRing3;
       phi0=BSDphi0Ring3;
       break;
     }
   phibase=phi-phi0-Twist*(z-BSDz0)/BSDscintHeight;
   if (phibase<=0.0) phibase+=360.0;
   if (phibase>360.0) phibase-=360.0;
   ScintNum=ceil((phibase*(float)NumScints)/360.0);
   number=ScintNum-1;
   energy[channel[number+chdim*key]]+=sorthits[4+nhdim*i]; 
   /* accumulate the energy */
   if ((tstart[channel[number+chdim*key]]+TDC_PULSE_RESOLUTION)<sorthits[6+nhdim*i]) { 
    /*start accumulating hits for this channel*/
    tstart[channel[number+key*chdim]]=sorthits[6+nhdim*i];
    Epulse[channel[number+key*chdim]]=sorthits[4+nhdim*i];
   }
   else {
    Epulse[channel[number+key*chdim]]+=sorthits[4+nhdim*i];
   }
   if (Epulse[channel[number+key*chdim]]>BSD_TDC_THRESHOLD) {
    /* front edge of the pulse */
    timelist[channel[number+chdim*key]+listdim*times[channel[number+chdim*key]]]=
                                         tstart[channel[number+chdim*key]];
    /* for the time-list */
    times[channel[number+chdim*key]]++; /* increase the times of hits saved */
    /* reset until next distinguishable pulse */
    Epulse[channel[number+chdim*key]]=-FLT_MAX;
   }
  }
 }
}

