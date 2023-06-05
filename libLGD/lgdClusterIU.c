/* lgdClusterIU.c 4.3 
   27 Aug 1996 22:29:13
   /auto/home/lemond/e852/source/Cluster/SCCS/s.lgdClusterIU.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <lgdGeom.h>
#include <lgdCluster.h>

static const char sccsid[] = "@(#)"__FILE__"\t4.3\tCreated 27 Aug 1996 22:29:13, \tcompiled "__DATE__;

static void levelOneSearch(lgd_hits_t *,const int16 *invertedIndex,lgd_clusters_t *,float, float *,lgd_hit_t **);
static int levelTwoSearch(lgd_hits_t *,lgd_clusters_t *,float,lgd_hit_t **,lgd_hits_t *);
static int getClusterPositions(lgd_clusters_t *,lgd_hit_t **);

#define SQR(x) ((x)*(x))
#define REALLOC(ptr,n) ((ptr)=realloc((ptr),(n)*sizeof((ptr)[0])))

static int level = LGDCLUS_DEF_LEVEL;
static float seedEnergy1 = LGDCLUS_DEF_SEEDENERGY1;
static float seedEnergy2 = LGDCLUS_DEF_SEEDENERGY2;
static int positionMethod = LGDCLUS_DEF_POSITION;
static float logOffset = LGDCLUS_DEF_LOGOFFSET;
static int maxNumClusters = LGDCLUS_DEF_MAXNUMCLUSTERS;
static int maxNumShared = LGDCLUS_DEF_MAXNUMSHARED;
static float maxRAssociate = LGDCLUS_DEF_MAXRASSOCIATE;
static float singleBlockWidth = LGDCLUS_DEF_SINGLEBLOCKWIDTH;

static int nChannels = 0;   /* The number of channels in the LGD */
static int nRows = 0;       /* The number of rows in the LGD */
static int nCols = 0;       /* The number of columns in the LGD */
static int *used = NULL;    /* Array, element = LGDCLUS_TRUE if 
			       the block of the same index in
			       <hits> is used in a cluster */

void lgdSetupAdvanced(int userLevel, float userSeedEnergy1,
		      float userSeedEnergy2, int userPositionMethod,
		      float userLogOffset,int userMaxNumClusters,
		      int userMaxNumShared,float userMaxRAssociate,
		      float userSingleBlockWidth)

     /* An advanced setup for the IU LGD clusterizer. 

	level is the level at which the clusterizer will operate -
	      1: Find only first generation clusters.
	      2: Find "associated" and "shared" blocks. Associated
	         blocks are blocks that were not put in a cluster
		 on the fist pass, but are within maxRAssociate cm
		 of one and only one cluster. Shared blocks are blocks
		 that were not put into a cluster on the first pass, but
		 are within maxRAssociate of two clusters. The energy of
		 the block is split between the two clusters. 
	      3: Find second generation showers.
	seedEnergy1 is the energy required in a single block to be the
	     seed for a first generation shower.
	seedEnergy2 is the energy required in a single block to be the
	     seed for a second generation shower.
	positionMethod tells the clusterizer what cluster position
	    finding alogorithm to use legal values are:
	        LGDCLUS_LINEAR for linear weighting
		LGDCLUS_LOG for log weighting
	logOffset is only used in if positionMethod is set to LGDCLUS_LOG.
	    It sets an offset that is used in determining the weight
	    given to a block.
	maxNumClusters is the maximum number of cluters that are allowed.
	    The clusterizer stops looking after it finds all that it is
	    allowed.
	maxNumShared is the maximum number of blocks that a cluster can
	    share with other clusters - This is only used in a level one
	    one search. An unlimited number of shared blocks are allowed
	    in the level 2 search.
	maxRAssociate is the maximum distance between a block and a cluster
	    to allow the block to be "associated" with the cluster in a level
	    2 search.
	singleBlockWidth is the width given to a cluster that contains only one
	    block 
     */


{
  level = userLevel;
  seedEnergy1 = userSeedEnergy1;
  seedEnergy2 = userSeedEnergy2;
  positionMethod = userPositionMethod;
  logOffset = userLogOffset;
  maxNumClusters = userMaxNumClusters;
  maxNumShared = userMaxNumShared;
  maxRAssociate = userMaxRAssociate;
  singleBlockWidth = userSingleBlockWidth;
}

void lgdSetup(int userLevel, float userSeedEnergy1, float userSeedEnergy2,
	      int userPositionMethod,float userLogOffset)

     /* Basic lgd setup. See above for more details */

{
  level = userLevel;
  seedEnergy1 = userSeedEnergy1;
  seedEnergy2 = userSeedEnergy2;
  positionMethod = userPositionMethod;
  logOffset = userLogOffset;
}

static void buildInvertedIndex(int nblocks,int16 *invertedIndex,lgd_hits_t *hits)
{
  int index;

  memset(invertedIndex,0,sizeof(invertedIndex[0])*nblocks);

  for(index=0;index<hits->nhits;index++)
    {
      if (hits->hit[index].channel < nblocks)
	{
	  invertedIndex[hits->hit[index].channel] = index + 1;
	}
    }
}

int lgdClusterizeIU(lgd_hits_t *hits,lgd_hits_t **userClusterHits,
		     lgd_clusters_t **userClusters)

     /* Clusterize the hits given in <hits>. Returns cluster info
	int <userClusterHits> and <clusters>.  These structures are
	correctly sized */


{
  float fractionClusterized=0.0;  /* The fraction of event energy used in
				     clusters */

  static int firstTime = 1;
  static int maxNhits = 0;
  static lgd_hit_t **hitList=NULL;       /* A 2-D array, containing the hits for
					    each cluster. hitList[0][1] is the
					    second hit in cluster 0 */
  static lgd_hits_t *isolatedHits=NULL;  /* The hits that were found to be 
					    isolated in the second level search */
  static lgd_clusters_t *clusters=NULL;  /* Temp. storage for cluster information */
  int index;                      /* A loop index */
  int nHitsUsed = 0;              /* The total number of hits clusterized */

  int16 invertedIndex[3053];      /* inverted index to map channel numbers
				     into entries in the hits array */

#ifdef UNDEF
  int16 invertedIndex1[3053];      /* inverted index to map channel numbers
				      into entries in the hits array */
#endif

  /* Get geometry info about the LGD */

  lgdGetGeomParam(&nChannels,&nRows,&nCols);

  /* Allocate memory */

  if ((hits->nhits > maxNhits)||firstTime)
    {
      firstTime = 0;

      maxNhits = hits->nhits;

      if(hitList){
	for(index=0;index<maxNumClusters;index++)
	  if(hitList[index])
	    free(hitList[index]);
	free(hitList);
      }

      hitList = malloc(sizeof(lgd_hit_t) * maxNumClusters);
      for(index=0;index<maxNumClusters;index++)
	hitList[index] = malloc(sizeof(lgd_hit_t)*hits->nhits);
  
      if(used)
	free(used);

      used = malloc(nChannels*sizeof(int));

      if(isolatedHits)
	free(isolatedHits);  

      isolatedHits = malloc((hits->nhits)*sizeof(lgd_hit_t) + sizeof(int));
      isolatedHits->nhits = 0;

      if(clusters)
	free(clusters);

      clusters = malloc(sizeof(lgd_clusters_t)*maxNumClusters+sizeof(int));
      clusters->nClusters = 0;
    }

  clusters->nClusters = 0;
  isolatedHits->nhits = 0;

  /* construct the inverted list of hits */

  buildInvertedIndex(3053,invertedIndex,hits);

  switch(level){
  case 1:
    levelOneSearch(hits,invertedIndex,clusters,seedEnergy1,&fractionClusterized,hitList);
    getClusterPositions(clusters,hitList);
    break;
  case 2:
    levelOneSearch(hits,invertedIndex,clusters,seedEnergy1,&fractionClusterized,hitList);
    getClusterPositions(clusters,hitList);
    levelTwoSearch(hits,clusters,fractionClusterized,hitList,isolatedHits);
    getClusterPositions(clusters,hitList);
    break;
  case 3:
    levelOneSearch(hits,invertedIndex,clusters,seedEnergy1,&fractionClusterized,hitList);
    getClusterPositions(clusters,hitList);
    levelTwoSearch(hits,clusters,fractionClusterized,hitList,isolatedHits);
    
    /* Search for clusters using the isolated hits list and seedEnergy2 */

#ifdef UNDEF
    buildInvertedIndex(3053,invertedIndex1,isolatedHits);
#endif
    levelOneSearch(isolatedHits,NULL /*invertedIndex1*/,clusters,seedEnergy2,&fractionClusterized,hitList);
    getClusterPositions(clusters,hitList);
    break;
  }
  
  /* Now allocate and fill userClusterHits */

  for(index=0;index<clusters->nClusters;index++)
    nHitsUsed += clusters->cluster[index].nBlocks;

  if(*userClusterHits)
    free(*userClusterHits);

  *userClusterHits = malloc(sizeof(lgd_hit_t)*nHitsUsed+sizeof(int));
  (*userClusterHits)->nhits = 0;

  for (index = 0; index < clusters->nClusters; index++)
    {
      memcpy(&(*userClusterHits)->hit[(*userClusterHits)->nhits],&hitList[index][0],
	     clusters->cluster[index].nBlocks*sizeof(lgd_hit_t));
      clusters->cluster[index].firstBlock = (*userClusterHits)->nhits;
      (*userClusterHits)->nhits += clusters->cluster[index].nBlocks;
    }
  
  /* Allocate and fill <userClusters */

  if(*userClusters)
    free(*userClusters);

  *userClusters = malloc(sizeof(lgd_cluster_t)*clusters->nClusters + sizeof(int));

  memcpy(*userClusters,clusters,sizeof(lgd_cluster_t)*clusters->nClusters + sizeof(int));

  /* Free some stuff up */

#ifdef UNDEF
  if(hitList){
    for(index=0;index<maxNumClusters;index++)
      if(hitList[index])
	free(hitList[index]);
    free(hitList);
  }
  hitList = NULL;

  if(used)
    free(used);
  used = NULL;

  if(isolatedHits)
    free(isolatedHits);
  isolatedHits = NULL;

  if(clusters)
    free(clusters);  
  clusters = NULL;
#endif

  return(1);
}


static void levelOneSearch(lgd_hits_t *hits,
		    const int16 *invertedIndex,
		    lgd_clusters_t *clusters,
		    float seedEnergy,float *fractionClusterized,
		    lgd_hit_t **hitList)

     /* Do a level 1 search for clusters in the list of hits <hits>.
	Use the seed energy specified in <seedEnergy>.  The cluster
	info. is put into <clusters>. The hits for each cluster are stored
	int hitList. <fractionClusterized> is set to the fraction of total
	event energy put into clusters */


{
  int index;                        /* A loop index */
  float eMax;                       /* The max. energy found in hits */
  int eMaxIndex;                    /* The index of eMax in hits */
  int neighbor[8];                  /* The neighbors of the max. energy hit */
  int nNeighbors;                   /* The number of neighbors */
  int neighborNum;                  /* An index for a loop over neighbors */

  static int clusterForHitAlloc = 0;
  static int *clusterForHit = NULL; /* An array, telling what cluster
				       eaxh hit is in, -1 -> not used yet */
  int nShared;                      /* The # of hits the current cluster
				       shares with other clusters */
  static int sharedWithAlloc = 0;
  static int *sharedWith = NULL;    /* The clusters the current cluster
				       shares the hit with */
  static int sharedHitAlloc = 0;
  static int *sharedHit = NULL;     /* The index in hitList of the shared
				       hits */
  int sharedHitNum;                 /* An index for a loop over shared hits */
  float totalEnergyClusterized;     /* The total energy in all clusters */
  static float totalEventEnergy;    /* The total energy in the event */
  float totalEnergyInCluster;       /* The total energy in the current 
				       cluster */
  int curClus;                      /* The current clusters index */
  int hitNum;                       /* An index, used in loops - either
				       the index in hits, or the index
				       in hitList */
  int curClusHitNum;                /* An index for a loop over the hits
				       in the current cluster */

  /* Allocate memory */

  if(sharedWithAlloc < maxNumShared)
    {
      sharedWithAlloc = maxNumShared;
      REALLOC(sharedWith,sharedWithAlloc);
    }

  if(sharedHitAlloc < maxNumShared)
    {
      sharedHitAlloc = maxNumShared;
      REALLOC(sharedHit,sharedHitAlloc);
    }

  if(clusterForHitAlloc < nChannels)
    {
      clusterForHitAlloc = nChannels;
      REALLOC(clusterForHit,clusterForHitAlloc);
    }

  if(clusters->nClusters == 0){
    
    /* Fist time through - find total energy in event */

    totalEventEnergy = 0.0;
    for(index=0;index<hits->nhits;index++)
      totalEventEnergy += hits->hit[index].energy;
  }

  for(index=0;index<hits->nhits;index++){
    
    /* mark all hits as unused */

    used[index] = LGDCLUS_FALSE;
    clusterForHit[index] = -1;
  }

  for(curClus=clusters->nClusters;curClus<maxNumClusters;curClus++){

    /* Find the block with the most energy that has not been used yet */

    eMax = -1.0;
    for(index=0;index<hits->nhits;index++)
      if (!used[index])
	if(hits->hit[index].energy > eMax)
	  {
	    eMax = hits->hit[index].energy;
	    eMaxIndex = index;
	  }
    
    if(eMax > seedEnergy){
     
      
      /* This block is a cluster seed -  put it into the structure */

      nShared = 0;

      memset(&clusters->cluster[clusters->nClusters],0, sizeof(clusters->cluster[clusters->nClusters]));

      clusters->nClusters++;      

      hitList[curClus][clusters->cluster[curClus].nBlocks++] = hits->hit[eMaxIndex];

      used[eMaxIndex] = LGDCLUS_TRUE;
      clusterForHit[eMaxIndex] = curClus;
      totalEnergyInCluster = hits->hit[eMaxIndex].energy;

      /* Get it's neighbors */      

      lgdGeomGetNeighbors(hits->hit[eMaxIndex].channel,neighbor,&nNeighbors);

      /* Check to see if any of the neighbors contain energy */

      for (neighborNum = 0; neighborNum < nNeighbors; neighborNum++){
	hitNum = (-1);

	if (invertedIndex != NULL)
	  {
	    hitNum = invertedIndex[neighbor[neighborNum]] - 1;
	  }
	else  /* no inverted index.... search through all the hits */
	  {
	    int i;

	    for (i = 0; i < hits->nhits; i++)
	      if(neighbor[neighborNum] == hits->hit[i].channel)
		{
		  hitNum = i;
		  break;
		}
	  }

	if (hitNum >= 0) { {

	    /* Neighbor has energy in it */

	    if(used[hitNum]){

	      /* It is already used in another cluster 
		 Put it into the current cluster, mark it as shared */

	      if(nShared < maxNumShared){
		hitList[curClus][clusters->cluster[curClus].nBlocks++] = hits->hit[hitNum];
		totalEnergyInCluster += hits->hit[hitNum].energy;
		sharedWith[nShared] = clusterForHit[hitNum];
		sharedHit[nShared] = hitNum;
		nShared++;
	      }
	    }
	    else{

	      /* It is not used in another cluster */

	      hitList[curClus][clusters->cluster[curClus].nBlocks++] = hits->hit[hitNum];
	      totalEnergyInCluster += hits->hit[hitNum].energy;
	      used[hitNum] = LGDCLUS_TRUE;
	      clusterForHit[hitNum] = curClus;
	    }
	  }
	}
      }
      clusters->cluster[curClus].energy = totalEnergyInCluster;
      
      if(nShared){

/*
	printf("nShared = %d\n",nShared);
*/
	
	/* Loop over shared hits */

	for(sharedHitNum = 0; sharedHitNum < nShared; sharedHitNum++) {

	  float fractInCurCluster;  /* The fraction of the energy of the
				       shared block that will go into the
				       current cluster */
	  float fractInOldCluster;  /* The fraction of the energy of the
				       shared block that will go into the
				       other cluster it is shared with */


	  fractInCurCluster = totalEnergyInCluster
	    /(totalEnergyInCluster + clusters->cluster[sharedWith[sharedHitNum]].energy);
	  fractInOldCluster = 1.0 - fractInCurCluster;

	  /* Loop over blocks in cluster that the hit is shared with */
	  
	  for(hitNum=0 ; 
	      hitNum < clusters->cluster[sharedWith[sharedHitNum]].nBlocks;
	      hitNum++){
	    if(hits->hit[sharedHit[sharedHitNum]].channel == 
	       hitList[sharedWith[sharedHitNum]][hitNum].channel){
	      
	      /* hitNum is the hit in the other cluster */
	      
	      for(curClusHitNum = 0;
		  curClusHitNum < clusters->cluster[curClus].nBlocks;
		  curClusHitNum++){
		if(hits->hit[sharedHit[sharedHitNum]].channel == 
		   hitList[curClus][curClusHitNum].channel){

		  /* curClusHitNum is the hit in this cluster 
		     share the energy between the two clusters */
		  
		  hitList[sharedWith[sharedHitNum]][hitNum].energy = fractInOldCluster * hitList[sharedWith[sharedHitNum]][hitNum].energy;
		  hitList[curClus][curClusHitNum].energy = fractInCurCluster* hitList[curClus][curClusHitNum].energy;
		}
	      }
	    }
	  }
	}

	/* Re-sum the energies in the clusters */

	clusters->cluster[curClus].energy = 0.0;
	for(hitNum = 0;
	    hitNum < clusters->cluster[curClus].nBlocks;
	    hitNum++)
	  clusters->cluster[curClus].energy += hitList[curClus][hitNum].energy;

	for (sharedHitNum = 0; sharedHitNum < nShared; sharedHitNum++) {
	  if(sharedHitNum > 0){
	    
	    /* Don't bother to resum the energies again if the last
	       shared clusters was the same as this one */

	    if(sharedWith[sharedHitNum-1] != sharedWith[sharedHitNum]){

	      /* Re-sum energy */

	      clusters->cluster[sharedWith[sharedHitNum]].energy = 0.0;
	      for(hitNum= 0; 
		  hitNum < clusters->cluster[sharedWith[sharedHitNum]].nBlocks;
		  hitNum++)
		clusters->cluster[sharedWith[sharedHitNum]].energy += 
		  hitList[sharedWith[sharedHitNum]][hitNum].energy;
	    }
	  }
	  else{

	    /* Re-sum energy */

	    clusters->cluster[sharedWith[sharedHitNum]].energy = 0.0;
	    for(hitNum=0; 
		hitNum < clusters->cluster[sharedWith[sharedHitNum]].nBlocks;
		hitNum++)
	      clusters->cluster[sharedWith[sharedHitNum]].energy +=
		hitList[sharedWith[sharedHitNum]][hitNum].energy;
	  }
	}
      }
    }
    else{
      
      /* No more clusters, quit now */

      totalEnergyClusterized = 0.0;
      for(index=0;index<clusters->nClusters;index++)
	totalEnergyClusterized += clusters->cluster[index].energy;
      *fractionClusterized = totalEnergyClusterized/totalEventEnergy;
      return;
      
    }
  }
  
  /* Have found all the clusters that are allowed */

  totalEnergyClusterized = 0.0;

  for (index = 0; index < clusters->nClusters; index++)
    totalEnergyClusterized += clusters->cluster[index].energy;

  (*fractionClusterized) = totalEnergyClusterized/totalEventEnergy;
  return;
}


static int levelTwoSearch(lgd_hits_t *hits,
		   lgd_clusters_t *clusters,
		   float fractionClusterized,lgd_hit_t **hitList,
		   lgd_hits_t *isolatedHits)

     /* Do a level 2 search of <hits> for "shared" and "associated" hits.
	Fill <isolatedHits> with the hits that remain unused */


{
  int hitNum;             /* An index, used in loop over hits */
  lgdGeom_t item;         /* The geometry info of the current hit */
  float rMin1,rMin2;      /* The square of the two shortest distances */
	
  float dist;             /* A temp variable, the square of the distance
			     between the current hit and the current
			     cluster */
  int cluster1,cluster2;  /* The indexes of the two closest clusters */
  int clusterNum;         /* An index for a loop over clusters */
  float fraction;         /* The fraction of the energy in the block that
			     will go to cluster1 */


  if(clusters->nClusters == 0 || fractionClusterized > 1.0)   return 0;

  for(hitNum=0;hitNum<hits->nhits;hitNum++){
    if(!used[hitNum]){
      
      /* Energy that is not in a cluster */

      lgdLocalCoord(hits->hit[hitNum].channel,&item.space);

      rMin1 = 999999.0;
      rMin2 = 999999.0;
      for(clusterNum=0;clusterNum<clusters->nClusters;clusterNum++){
	dist = SQR(item.space.x - clusters->cluster[clusterNum].space.x) +
	       SQR(item.space.y - clusters->cluster[clusterNum].space.y);
	if (dist < rMin1) {
          rMin2 = rMin1;
          cluster2 = cluster1;
	  rMin1 = dist;
	  cluster1 = clusterNum;
	}
	else if (dist < rMin2) {
	  rMin2 = dist;
	  cluster2 = clusterNum;
	}
      }

      if( rMin1 <= SQR(maxRAssociate) && rMin2 > SQR(maxRAssociate)){
	
	/* This block is near only one cluster - add it to the cluster */
	/* This is a associated block */

	used[hitNum] = LGDCLUS_TRUE;
	hitList[cluster1][clusters->cluster[cluster1].nBlocks++] = hits->hit[hitNum];
	clusters->cluster[cluster1].energy += hitList[cluster1][clusters->cluster[cluster1].nBlocks -1].energy;
      }
      else{
	if(rMin2 <= SQR(maxRAssociate)){
	  
	  /* Block is shared */

	  used[hitNum] = LGDCLUS_TRUE;
	  fraction = clusters->cluster[cluster1].energy/
	    (clusters->cluster[cluster1].energy+clusters->cluster[cluster2].energy);
	  hitList[cluster1][clusters->cluster[cluster1].nBlocks++] = hits->hit[hitNum];
	  hitList[cluster2][clusters->cluster[cluster2].nBlocks++] = hits->hit[hitNum];
	  hitList[cluster1][clusters->cluster[cluster1].nBlocks -1].energy *= fraction;
	  hitList[cluster2][clusters->cluster[cluster2].nBlocks -1].energy *= (1.0 - fraction);
	  clusters->cluster[cluster1].energy += hitList[cluster1][clusters->cluster[cluster1].nBlocks -1].energy;
	  clusters->cluster[cluster2].energy += hitList[cluster2][clusters->cluster[cluster2].nBlocks -1].energy;
	}
	else{

	  /* Block is isolated */
	  /* Put it into new "isolated" hit list */

	  isolatedHits->hit[isolatedHits->nhits++] = hits->hit[hitNum];
	}
      }
    }
  }
  return 0;
}


static int getClusterPositions(lgd_clusters_t *clusters,lgd_hit_t **hitList)

     /* Find the X,Y,and Z coord of the clusters in LGD space */

{
  int clusterNum;      /* An index for a loop over clusters */
  int blockNum;        /* An index for a loop over blocks in the cluster */
  lgdGeom_t  block;    /* The geometry info for the current block */
  float currentOffset; /* The current log offset for the block */
  lgdGeom_t centerBlock; /* What it sounds like... */

  /* For LOG2 */

  float centerWeight;
  float neighborMaxWeight;
  float logFraction;
  

  /* Stuff for "Bugge" method */

  float maxNeighborEx;
  float signMaxNeighborEx;
  float maxNeighborEy;
  float signMaxNeighborEy;
  float correctionThreshold = BUGGE_MAGIC_THRESHOLD;
  float sigmaA=logOffset;


  switch(positionMethod){
  case LGDCLUS_LINEAR:

    /* Use "linear" weighting */

    for(clusterNum=0;clusterNum < clusters->nClusters;clusterNum++){
      clusters->cluster[clusterNum].space.x = 0.0;
      clusters->cluster[clusterNum].space.y = 0.0;
      clusters->cluster[clusterNum].space.z = 0.0;
      for(blockNum=0;blockNum<clusters->cluster[clusterNum].nBlocks;blockNum++){
	lgdLocalCoord(hitList[clusterNum][blockNum].channel,&block.space);
	clusters->cluster[clusterNum].space.x += block.space.x*hitList[clusterNum][blockNum].energy;
	clusters->cluster[clusterNum].space.y += block.space.y*hitList[clusterNum][blockNum].energy;
      }
      clusters->cluster[clusterNum].space.x /= clusters->cluster[clusterNum].energy;
      clusters->cluster[clusterNum].space.y /= clusters->cluster[clusterNum].energy;
    }
    break;
  case LGDCLUS_LOG:

    /* Use "log" weighting */

    for(clusterNum=0;clusterNum < clusters->nClusters;clusterNum++){
      float weight = 0.0;
      float weightSum = 0.0;
      clusters->cluster[clusterNum].space.x = 0.0;
      clusters->cluster[clusterNum].space.y = 0.0;
      clusters->cluster[clusterNum].space.z = 0.0;
      currentOffset = clusters->cluster[clusterNum].energy/3.0 + 2.6;
      currentOffset = (currentOffset > logOffset)?logOffset:currentOffset;
      for(blockNum=0;blockNum<clusters->cluster[clusterNum].nBlocks;blockNum++){
	lgdLocalCoord(hitList[clusterNum][blockNum].channel,&block.space);
	/* Find the weight - offset determines minimun energy of block
	   to be used in weighting, since negative weights are thrown out */

	weight = currentOffset + log(hitList[clusterNum][blockNum].energy/clusters->cluster[clusterNum].energy);
	if(weight > 0.0){
	  clusters->cluster[clusterNum].space.x += block.space.x*weight;
	  clusters->cluster[clusterNum].space.y += block.space.y*weight;
	  weightSum += weight;
	}
      }
      clusters->cluster[clusterNum].space.x /= weightSum;
      clusters->cluster[clusterNum].space.y /= weightSum;
    }
    break;
  case LGDCLUS_LOG2:

    /* Use "log" weighting */

    for(clusterNum=0;clusterNum < clusters->nClusters;clusterNum++){
      float weight = 0.0;
      float weightSum = 0.0;
      clusters->cluster[clusterNum].space.x = 0.0;
      clusters->cluster[clusterNum].space.y = 0.0;
      clusters->cluster[clusterNum].space.z = 0.0;
      neighborMaxWeight=0.0;
      logFraction = exp(-0.23*(clusters->cluster[clusterNum].energy));
      currentOffset = log(clusters->cluster[clusterNum].energy)/7.0 + 3.7;
      for(blockNum=0;blockNum<clusters->cluster[clusterNum].nBlocks;blockNum++){
	lgdLocalCoord(hitList[clusterNum][blockNum].channel,&block.space);

	/* Find the weight - offset determines minimun energy of block
	   to be used in weighting, since negative weights are thrown out */

	weight = currentOffset + log(hitList[clusterNum][blockNum].energy/clusters->cluster[clusterNum].energy);
	if(blockNum==0){
	  centerWeight = weight;
	  lgdLocalCoord(hitList[clusterNum][blockNum].channel,&centerBlock.space);
	}
	else{
	  if(neighborMaxWeight<weight)
	    neighborMaxWeight= weight;
	}
	if(weight > 0.0){
	  clusters->cluster[clusterNum].space.x += block.space.x*weight;
	  clusters->cluster[clusterNum].space.y += block.space.y*weight;
	  weightSum += weight;
	}
      }

      /* now patch up the center block's weight if it's got a neighbor in the cluster that had positive weight*/

      if(neighborMaxWeight > 0){
	clusters->cluster[clusterNum].space.x+=(logFraction-1)*(centerWeight-neighborMaxWeight)*centerBlock.space.x;
	clusters->cluster[clusterNum].space.y+=(logFraction-1)*(centerWeight-neighborMaxWeight)*centerBlock.space.y;
	weightSum += (logFraction-1)*(centerWeight-neighborMaxWeight);
      }

      clusters->cluster[clusterNum].space.x /= weightSum;
      clusters->cluster[clusterNum].space.y /= weightSum;
    }
    break;
  case LGDCLUS_BUGGE:

    /* Use "Bugge" weighting */


    for(clusterNum=0;clusterNum < clusters->nClusters;clusterNum++){
      clusters->cluster[clusterNum].space.x = 0.0;
      clusters->cluster[clusterNum].space.y = 0.0;
      clusters->cluster[clusterNum].space.z = 0.0;

      maxNeighborEx=0.0;
      signMaxNeighborEx=0.0;
      maxNeighborEy=0.0;
      signMaxNeighborEy=0.0;

      for(blockNum=0;blockNum<clusters->cluster[clusterNum].nBlocks;blockNum++){
        lgdLocalCoord(hitList[clusterNum][blockNum].channel,&block.space);


        /* some book keeping stuff for later use */
        if(blockNum==0){
          lgdLocalCoord(hitList[clusterNum][blockNum].channel,&centerBlock.space);
        }
        else{
          if(block.space.x==centerBlock.space.x){
            if((block.space.y==(centerBlock.space.y+4.0))&&(hitList[clusterNum][blockNum].energy>maxNeighborEy)){
              maxNeighborEy=hitList[clusterNum][blockNum].energy;
              signMaxNeighborEy=+1.0;
            }
            else{
              if((block.space.y==(centerBlock.space.y-4.0))&&(hitList[clusterNum][blockNum].energy>maxNeighborEy)){
                maxNeighborEy=hitList[clusterNum][blockNum].energy;
                signMaxNeighborEy=-1.0;
              }
            }
          }
          
          if(block.space.y==centerBlock.space.y){
            if((block.space.x==(centerBlock.space.x+4.0))&&(hitList[clusterNum][blockNum].energy>maxNeighborEx)){
              maxNeighborEx=hitList[clusterNum][blockNum].energy;
              signMaxNeighborEx=+1.0;
            }
            else{
              if((block.space.x==(centerBlock.space.x-4.0))&&(hitList[clusterNum][blockNum].energy>maxNeighborEx)){
                maxNeighborEx=hitList[clusterNum][blockNum].energy;
                signMaxNeighborEx=-1.0;
              }
            }
          }
        }
        
        clusters->cluster[clusterNum].space.x += block.space.x*hitList[clusterNum][blockNum].energy;
        clusters->cluster[clusterNum].space.y += block.space.y*hitList[clusterNum][blockNum].energy;
      }
      clusters->cluster[clusterNum].space.x /= clusters->cluster[clusterNum].energy;
      clusters->cluster[clusterNum].space.y /= clusters->cluster[clusterNum].energy;

     
      /* okay now let's bastardize it with some non-linear corrections if there was neighbor we can use*/
      /* first do x */
      if(maxNeighborEx>correctionThreshold){
        clusters->cluster[clusterNum].space.x = centerBlock.space.x + signMaxNeighborEx*2.0 - 
          signMaxNeighborEx*sigmaA*(log(0.5*(hitList[clusterNum][0].energy/maxNeighborEx +1)));
      }
      /* now do y */
      if(maxNeighborEy>correctionThreshold){
        clusters->cluster[clusterNum].space.y = centerBlock.space.y + signMaxNeighborEy*2.0 - 
          signMaxNeighborEy*sigmaA*(log(0.5*(hitList[clusterNum][0].energy/maxNeighborEy +1)));
      }  




    }
    break;
  }
  
  /* Calculate energy weighted second moment, and use it as the width of the
     cluster */

  for(clusterNum=0;clusterNum<clusters->nClusters;clusterNum++){
    clusters->cluster[clusterNum].width = 0.0;
    if(clusters->cluster[clusterNum].nBlocks == 1)

      /* Use single block cluster width */

      clusters->cluster[clusterNum].width = singleBlockWidth;
    else{
      for(blockNum=0;blockNum<clusters->cluster[clusterNum].nBlocks;blockNum++){
	lgdLocalCoord(hitList[clusterNum][blockNum].channel,&block.space);
	clusters->cluster[clusterNum].width += hitList[clusterNum][blockNum].energy*
	  (SQR(block.space.x - clusters->cluster[clusterNum].space.x) +
	   SQR(block.space.y - clusters->cluster[clusterNum].space.y));
      }
      clusters->cluster[clusterNum].width = sqrt(clusters->cluster[clusterNum].width/
						 clusters->cluster[clusterNum].energy);
    }
  }    
  return 0;
}
