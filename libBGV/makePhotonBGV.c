/* returns the 4-vector of a photon that created cluster */

// #define BGV_NONLINEARITY_CORRECTION 1  /* disabled by default -rtj */

#include <math.h>

#include <clib.h>
#include <itypes.h>
#include <iitypes.h>
#include <disData.h>
#include <libBGV.h>
#include "bgvGeom.h"

vector4_t makePhotonBGV(const bgv_cluster_t *cluster)
{
  float cosTheta,sinTheta;
  float cosPhi,sinPhi;
  float px,py,pz,pmag;
  float energy,theta;
  vector4_t p;

#ifdef VERBOSE
  if (cluster->energy > 1e9) 
    fprintf(stderr,"Warning: very high energy in makePhotonBGV (E > 1 billion GeV)\n");
  if (cluster->energy == 0) {
    fprintf(stderr,"warning: energy of cluster is 0\n");
    fprintf(stderr,"  z=%f, phi=%f\n",cluster->z,cluster->phi);
  }
  if ((cluster->z < -100)||(cluster->z > 200)) {
    fprintf(stderr,"warning: photon not in BGV in function makePhotonBGV, z=%f\n",
	   cluster->z);
  }
#endif
  energy=cluster->energy;
  theta=atan2(MEAN_RADIUS_BGV,cluster->z);

#ifdef BGV_NONLINEARITY_CORRECTION

/* This new nonlinearity correction for the barrel was based
 * on a study of single-photon generation Monte Carlo events.
 * Richard.T.Jones@uconn.edu, October 4, 2005
 */

  {
    float par0[] = {0.24765,-0.75997,0.46902,-0.09652};
    float par1[] = {0.89199,0.036188,-0.30577,0.12660};
    float par2[] = {-0.48,0.48};

    float Ecut = 3.0; // GeV
    float Tcut = 2.0; // radians
    float resolution = 0.001;
    float Eorig = energy;
    float Torig = theta;
    int niter;
    for (niter=0; 1; niter++) {
      float p0 = (((par0[3])*Torig+par0[2])*Torig+par0[1])*Torig+par0[0];
      float p1 = (((par1[3])*Torig+par1[2])*Torig+par1[1])*Torig+par1[0];
      float p2 = (Eorig-par2[0])/par2[1];
      float Eobs = Eorig*(p0*sqrt(p2-1)+p1);
      float Tobs = Torig;
      if (fabs(Eobs-energy)+fabs(Tobs-theta) < resolution) {
        energy = Eorig;
        theta = Torig;
        break;
      }
      else if (Eorig > Ecut || theta > Tcut) {
        energy = 0;
        break;
      }
      else if (niter > 1000) {
        fprintf(stderr,"makePhotonBGV - maximum iteration count exceeded,");
        fprintf(stderr," giving up!\n");
        energy = 0;
        break;
      }
      Eorig *= energy/Eobs;
      Torig *= theta/Tobs;
    }
    if ((theta != theta) || (energy != energy)) {
      fprintf(stderr,"makePhotonBGV - iterative solution failed!\n");
    }
  }

#endif

  cosTheta=cos(theta);
  sinTheta=sin(theta);
  cosPhi=cos(cluster->phi);
  sinPhi=sin(cluster->phi);
  if (energy > 0) {
    pmag=energy;
  }
  else {
    pmag=cluster->energy;
  }
  p.space.z=pmag*cosTheta;
  p.space.x=pmag*sinTheta*cosPhi;
  p.space.y=pmag*sinTheta*sinPhi;
  p.t=energy;
  return(p);
}
