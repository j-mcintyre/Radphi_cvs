/*     makePhotonPatch.c
 *
 *     D.S. Armstrong    March 9 2000
 *     - original version from C. Steffen
 *
 *     R.T. Jones	 February 2 2001
 *     - replaced old depth correction with theta-dependent correction
 *       that solves for corrected energy and angle at the same time.
 *       For more information see TechNote 2001-201.
 *
 *     R.T. Jones	 Spetember 28 2005
 *     - update the theta-dependent depth and nonlinearity correction
 *     - for clusters with energy,theta combinations that are outside the
 *       range of allowable reconstruction resolution (inside and outside
 *       boundaries) replace the corrected shower energy with 0.
 */

//#define SIMPLE_DEPTH_CORRECTION -nominal depth formula, no energy correction
//#define OLD_DEPTH_CORRECTION    -original depth correction

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>  
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <ntypes.h>
#include <dataIO.h>
#include <itypes.h>
#include <iitypes.h>
#include <mctypes.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <lgdGeom.h>
#include <disData.h>
#include <pedestals.h>
#include <calibration.h>
#include <disIO.h>
#include <eventType.h>
#include <triggerType.h>
#include <makeHits.h>

vector4_t makePhotonPatch(lgd_cluster_t *cluster,vector3_t vertex){

#define LGD_RADIATION_LENGTH 3.1
#define LGD_SHOWER_OFFSET 3.0
#define E_CRITICAL 0.01455
#define RAD_DEG 57.29578

#define SQR(X) ((X)*(X))

  vector3_t photon_VC; /* photon in vertex coordinates */
  vector4_t answer;
#ifdef OLD_DEPTH_CORRECTION
  float radius,baseRadius,newRadius,dr;
  float zMax;
#endif
  int niter = 0;
  float energy;
  float magnitude;
  float scale;

/* The local coordinate system for the LGD is situated with the origin
 * centered on the entrance face to the lead glass array, with the axes
 * parallel to the global coordinate axes.
 */

  lgdLocalToGlobalCoord(cluster->space,&photon_VC);
  energy = cluster->energy;
  photon_VC.x -= vertex.x;
  photon_VC.y -= vertex.y;
  photon_VC.z -= vertex.z;

#if defined SIMPLE_DEPTH_CORRECTION

  photon_VC.z += LGD_RADIATION_LENGTH*(2.0+log(energy/E_CRITICAL));

#elif defined OLD_DEPTH_CORRECTION

/* The following code was formerly used to estimate the z coordinate of
 * the shower centroid in the LGD, under the assumption that the wall
 * is infinite in x and y.  It has been replaced with a treatment that
 * depends on polar angle -- for more information see TechNote 2001-201.
 */

  /* setting initial conditions */
  radius=sqrt(SQR(photon_VC.x)+SQR(photon_VC.y));
  baseRadius = radius;
  newRadius = radius;
 
  /* Get the depth of the shower max */
  zMax = LGD_RADIATION_LENGTH*
         (LGD_SHOWER_OFFSET + log(cluster->energy/E_CRITICAL));

  dr = 99.0;
  while(fabs(dr) > 0.001){
    newRadius = baseRadius/(1+(zMax/sqrt(SQR(photon_VC.z) + SQR(radius))));
    dr = newRadius-radius;
    radius = newRadius;
    niter++;
    if (niter > 20){
      printf("Help!  I'm stuck inside a while loop with %d iterations!",niter);
      break;
    }
  }
  /* printf("old: %f\t%f\t%f\n",photon_VC.x,photon_VC.y,photon_VC.z); */
  photon_VC.x *= (newRadius/baseRadius);
  photon_VC.y *= (newRadius/baseRadius);
  /* printf("new: %f\t%f\t%f\n",photon_VC.x,photon_VC.y,photon_VC.z); */

#else
Iterator_start:

/* Start of iterative loop where corrected energy,angle are solved in a
 * self-consistent way -- for more information see TechNote 2001-201 --rtj.
 *
 * -modified 27-Sep-2005 [rtj]
 *  This code was so messed up, with patches on top of patches, that I
 *  just cleaned it out and replaced it with a fresh treatment.  There is
 *  no longer any attempt to recover showers outside the range of validity
 *  of the energy/depth correction.  In that case, just return a nominal
 *  shower centroid and a corrected energy of zero.
 */

  {
    float zed = photon_VC.z+LGD_RADIATION_LENGTH*(2.0+log(energy/E_CRITICAL));
    float theta = atan2(sqrt(SQR(photon_VC.x)+SQR(photon_VC.y)),zed);
    float par0[] = {1.5479,-16.249,193.35,-1140.7,3472.3,-5383.1,3424.1};
    float par1[] = {-0.0018,-0.0136,-25.213,54.761};
    float par2[] = {1.0590,-0.60352,2.4527,-3.8784};
    float par3[] = {-0.00418,0.03098,-0.0533,-0.0554};

    float Ecut = 8.0; // GeV
    float resolution = 0.001;
    float Eorig = energy;
    float Torig = theta;
    for (niter=0; 1; niter++) {
      float p0 = ((((((par0[6])*Torig+par0[5])*Torig+par0[4])*Torig
                      +par0[3])*Torig+par0[2])*Torig+par0[1])*Torig+par0[0];
      float p1 = par1[0]+par1[1]*Torig+exp(par1[2]+par1[3]*Torig);
      float Eobs = Eorig/(p0+p1*Eorig);
      float p2 = ((par2[3]*Torig+par2[2])*Torig+par2[1])*Torig+par2[0];
      float p3 = ((par3[3]*Torig+par3[2])*Torig+par3[1])*Torig+par3[0];
      float Tobs = Torig*(p2+p3*Eorig);
      if (fabs(Eobs-energy)+fabs(Tobs-theta) < resolution) {
        energy = Eorig;
        theta = Torig;
        break;
      }
      else if (Eorig > Ecut) {
        energy = 0;
        break;
      }
      else if (niter > 1000) {
        fprintf(stderr,"makePhotonPatch - maximum iteration count exceeded,");
        fprintf(stderr," giving up!\n");
        energy = 0;
        break;
      }
      Eorig *= energy/Eobs;
      Torig *= theta/Tobs;
    }
    if ((theta != theta) || (energy != energy)) {
      fprintf(stderr,"makePhotonPatch - iterative solution failed!\n");
    }
    photon_VC.z = sqrt(SQR(photon_VC.x)+SQR(photon_VC.y))/tan(theta);
  }

#endif

  answer.space.x=photon_VC.x;
  answer.space.y=photon_VC.y;
  answer.space.z=photon_VC.z;
  answer.t = energy;

  /* answer.space now has the proper angles; now rescale it so that those 
     mean momenta */

  magnitude=sqrt(SQR(answer.space.x)+SQR(answer.space.y)+SQR(answer.space.z));
  scale=answer.t/magnitude;
  if (scale > 0) {
    answer.space.x *= scale;
    answer.space.y *= scale;
    answer.space.z *= scale;
  }
  return answer;
}
