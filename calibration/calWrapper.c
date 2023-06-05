#include <lgdCluster.h>
#include <lgdUtil.h>
#include <math.h>
#include <stdio.h>
#include <clib.h>
#include <lgdCal.h>


void ecal0_(int ncc[CAL_PRIM_SIZE],float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],float x[CAL_PRIM_SIZE],float u[CAL_PRIM_SIZE]);
void ecal_(float *efit,int nbclus[25][8],float ebclus[25][8],int mbclus[8],
	   int *index,
	   float cc[CAL_PRIM_SIZE],int ncc[CAL_PRIM_SIZE], float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],
	   float x[CAL_PRIM_SIZE],float u[CAL_PRIM_SIZE]);
void esolv_(float cc[CAL_PRIM_SIZE] ,float newcc[CAL_PRIM_SIZE], int ncc[CAL_PRIM_SIZE],
	    float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],float x[CAL_PRIM_SIZE],float u[CAL_PRIM_SIZE]);

void ecal0(int ncc[CAL_PRIM_SIZE],float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],float x[CAL_PRIM_SIZE],
	   float u[CAL_PRIM_SIZE])

{
  ecal0_(ncc,store,x,u);
}

void ecal(float efit,lgd_cluster_t *cluster,lgd_hits_t *hits,
	  float cc[CAL_PRIM_SIZE],int ncc[CAL_PRIM_SIZE],float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],
	  float x[CAL_PRIM_SIZE],float u[CAL_PRIM_SIZE])

{
  /* Now put the cluster info into Scotts "format" */

  int nbclus[25][8];
  float ebclus[25][8];
  int mbclus[8]={0,0,0,0,0,0,0,0};
  int index;

  if( (efit < 0) || (fabs(efit) < 0.00001) || (is_NaN(efit)))
    fprintf(stderr,"Bad fitted energy!\n");
  
  for(index=cluster->firstBlock;index<cluster->firstBlock+cluster->nBlocks;index++){
    
    /* fill nbclus - add one to channel to account for Scotts
       variation on the indexing system */

    nbclus[mbclus[0]][0] = hits->hit[index].channel+1;
    ebclus[mbclus[0]][0] = hits->hit[index].energy;
    mbclus[0]++;
  }
  
  /* Tell him it's cluster #1 */

  index = 1;
  ecal_(&efit,nbclus,ebclus,mbclus,&index,cc,ncc,store,x,u);
}

void esolv(float cc[CAL_PRIM_SIZE],float newcc[CAL_PRIM_SIZE],int ncc[CAL_PRIM_SIZE],
	   float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],float x[CAL_PRIM_SIZE],float u[CAL_PRIM_SIZE])

{
  esolv_(cc,newcc,ncc,store,x,u);
}
