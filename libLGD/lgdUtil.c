/*
  $Log: lgdUtil.c,v $
  Revision 1.2  1997/05/29 20:44:17  radphi
  Changes made by lfcrob@dustbunny
  Removed utility.h header

 * Revision 1.1  1997/05/24  00:57:17  radphi
 * Initial revision by lfcrob@jlabs2
 * Misc. LGD function
 *
  */

static const char rcsid[] = "$Id: lgdUtil.c,v 1.2 1997/05/29 20:44:17 radphi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>

#include <itypes.h>
#include <lgdCluster.h>
#include <lgdGeom.h>
#include <lgdUtil.h>
#include <clib.h>

#define sqr(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define DEMIN 0.00001
#define MAXCOUNT 10
#define E_CRITICAL 0.01455
#define ZMAX_COEFF 4.52

#define SQR(x) ((x)*(x))

static float depthScale;

vector4_t makePhoton(lgd_cluster_t *cluster, vector3_t vertex)
{
  return(makePhotonGen(MAKE_PHOTON_DEFAULT,cluster,NULL,vertex));
}

vector4_t makePhotonGen(int mode,lgd_cluster_t *cluster, lgd_hits_t *clusterHits, vector3_t vertex)
{
  lgd_cluster_t oldCluster;
  float eprime;
  float de;
  float scale;
  vector4_t p;
  int i;
  int depthMode=LGD_DEPTH_NONE;

  if(mode & MAKE_PHOTON_DEPTH){
    depthMode = LGD_DEPTH_1;
    if(!(mode & MAKE_PHOTON_USR_DEPTH))  
      setupClusterToPhoton(0.8);
  }  
  p = clusterToPhoton1(depthMode,cluster,vertex);
  return(p);
}

void makePhotons(lgd_clusters_t *clusters,lgd_hits_t *clusterHits,vector4_t *p,vector3_t vertex)
{
  int index;
  
  for(index=0;index<clusters->nClusters;index++){
    p[index]=makePhotonGen(MAKE_PHOTON_DEFAULT,&clusters->cluster[index],clusterHits,vertex);
  }
}

void setupClusterToPhoton(float scale)
{
  depthScale = scale;
}

vector4_t clusterToPhoton(lgd_cluster_t *cluster,vector3_t vertex)
{
  return clusterToPhoton1(LGD_DEPTH_1,cluster,vertex);
}


vector4_t clusterToPhoton1(int mode,lgd_cluster_t *cluster,vector3_t vertex)
{
  vector3_t globalVect,vertexVect;
  vector4_t answer;
  float mag;

  float radius,baseRadius,newRadius,dr;
  float zMax;

  /*Get the cluster position in global Coord */

  lgdLocalToGlobalCoord(cluster->space,&globalVect);

  /* Get the cluster position in "vertex" coord */

  v3diff(&globalVect,&vertex,&vertexVect);

  /* Get the cluster radius in vertex coord. */

  radius = sqrt(sqr(vertexVect.x)+sqr(vertexVect.y));

  baseRadius = radius;
  newRadius = radius;

  switch (mode)
    {
    default:
    case LGD_DEPTH_NONE:
      break;

    case LGD_DEPTH_1:
      { 
	int niter = 0;
	/* Get the depth of the shower max */
	
	zMax = depthScale*ZMAX_COEFF*log(cluster->energy/E_CRITICAL);
	
	dr = 99.0;
	while(fabs(dr) > 0.001){
	  newRadius = baseRadius/(1+(zMax/
				     sqrt(sqr(vertexVect.z) + sqr(radius))));
	  dr = newRadius-radius;
	  radius = newRadius;
	  niter ++;
	  if (niter > 20)
	    break;
	}
      }
      break;
    }
  
  /*Rescale the Target coord cluster pos, go to Global coord */

  vertexVect.y *= (newRadius/baseRadius);
  vertexVect.x *= (newRadius/baseRadius);

  answer.space = vertexVect;
  
  /* Get the gamma 3-vector mag */

  mag = v3mag(&answer.space);

  /* Scale the photon 4-vector */

  answer.t = cluster->energy;
  ScalerMult(&answer.space,cluster->energy/mag,&answer.space);
  return(answer);
}


float twoGammaMass(vector4_t *p1,vector4_t *p2)
{
  vector4_t *pArray[2];
  int nPart=2;
  
  pArray[0] = p1;
  pArray[1] = p2;
  return(EffMass(nPart,(const vector4_t **)pArray));  
}

float fitTwoGamma(vector4_t *p1,vector4_t *p2,float mass,
		  float sigma1,float sigma2,
		  vector4_t *pFit1, vector4_t *pFit2)

     /* Do constrained fit of two gammas to <mass>
	Return value is the chisq of the fit */

{
  double a,b,c,q;
  double de;
  int index;
  float chisq;
  int count;

  de = 99999.0;
  count=0;
  pFit1->t = p1->t;
  q = sqr(mass)/(2*(1-ctheta(&p1->space,&p2->space)));
  a = -p1->t;
  b = q*p2->t*sqr(sigma1/sigma2);
  c = -sqr(q*(sigma1/sigma2));
  while((fabs(de) > DEMIN) && (++count < MAXCOUNT)){
    de = (c+pFit1->t*(b+sqr(pFit1->t)*(a+pFit1->t)))/
      (b+sqr(pFit1->t)*(3*a+4*pFit1->t));
    pFit1->t -= de;
  }
  pFit2->t = q/pFit1->t;
  ScalerMult(&pFit1->space,pFit1->t/p1->t,&p1->space);
  ScalerMult(&pFit2->space,pFit2->t/p2->t,&p2->space);
  chisq = sqr((p1->t-pFit1->t)/sigma1)+sqr((p2->t-pFit2->t)/sigma2);
  return(chisq);
}
