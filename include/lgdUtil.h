/*
  
  $Log: lgdUtil.h,v $
  Revision 1.2  2000/03/13 21:56:02  radphi
  Changes made by armd@jlabs1
  change MAKE_PHOTON_USR_DEPTH so that default depthScale is not automatically overwritten

  Revision 1.1  1997/05/17 19:07:00  radphi
  Initial revision by lfcrob@jlabs2
  Misc LGD routines

*/


#ifndef LGDUTIL_H_INCLUDED
#define LGDUTIL_H_INCLUDED


#define MAKE_PHOTON_DEPTH 0x1
#define MAKE_PHOTON_USR_DEPTH 0x1
/*  changed dsa */
/* #define MAKE_PHOTON_USR_DEPTH 0x2 */
#define MAKE_PHOTON_DEFAULT MAKE_PHOTON_DEPTH

vector4_t makePhoton(lgd_cluster_t *cluster,vector3_t vertex);
vector4_t makePhotonGen(int mode,lgd_cluster_t *cluster,lgd_hits_t *clusterHits,vector3_t vertex);
void makePhotons(lgd_clusters_t *clusters,lgd_hits_t *clusterHits,vector4_t *p,vector3_t vertex);


void setupRescale2(float,float,float,float,float);


/* Old (IU-1994) LGD cluster correction code */

#define LGD_DEPTH_NONE    0
#define LGD_DEPTH_1      0x1

/*
 * clusterToPhoton1(...) returns the gamma four-vector corresponding
 * to an LGD cluster. Depending on the 'mode', depth correction
 * and correction for non-normal incidence might be used.
 *
 * clusterToPhoton(...) is identical to clusterToPhoton1(...)
 * with 'mode' equal to 'LGD_ClusterCorrection_IU'
 * ONLY FOR OLD CC's!!!!!! (Pre-95)
*/

void setupClusterToPhoton(float scale);    /* setup the Z-correction for the LGD_ClusterCorrection_IU94 mode */
vector4_t clusterToPhoton(lgd_cluster_t *cluster,vector3_t vertex);   /* use the LGD_ClusterCorrection_IU94 mode */
vector4_t clusterToPhoton1(int mode, lgd_cluster_t *cluster,vector3_t vertex);   /* use the LGD_ClusterCorrection_IU94 mode */


/*
 * perform the constrained fit of the two gamma effective mass
*/

float fitTwoGamma(vector4_t *p1,vector4_t *p2,float mass,
		  float sigma1,float sigma2,
		  vector4_t *pFit1, vector4_t *pFit2);

/*
 * calculate the effective mass of 'npart' four vectors
*/


#define SIGMA_THETA_DEF 0.00116

#define FIT_ENERGY_ONLY 1
#define FIT_POSITION_1 2


int fitTwoGammaFull(vector3_t vertex,vector4_t *p1,vector4_t *p2,float mass,
		    float sigma1,float sigma2,float sigmaTheta,
		    vector4_t *p1Fit, vector4_t *p2Fit,float *chisq);
int fitTwoGammaGen(int mode,
		   vector3_t vertex,vector4_t *p1,vector4_t *p2,float mass,
		   float sigma1,float sigma2,float sigmaTheta,
		   vector4_t *p1Fit, vector4_t *p2Fit,float *chisq);

#endif
