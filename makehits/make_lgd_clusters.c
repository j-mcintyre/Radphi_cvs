/*     make_lgd_clusters.c     */
/*     D.S. Armstrong    March 9 2000  */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <itypes.h> 
#include <eventType.h>
#include <ntypes.h> 

#include <disData.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <lgdGeom.h>
#include <pedestals.h>
#include <calibration.h> 
#include <makeHits.h>

#define MAX_PATH 256
//#define SPLITOFF_RECLAMATION true
#define CLUSTER_RECLAIM_THRESHOLD 6
//#define PRUNE_SINGLES true
//#define CLUSTER_CLEANUP_SEPARATION_CUT 10  
//#define SINGLE_BLOCK_CLUSTER_THRESHOLD .150
#define SQR(x) ((x)*(x))

/* variables global in this object file.  */
int lgd_cluster_cleanup_initialized=0;
float lgd_min_cluster_rho; 
  /* smallest cluster distance from LGD center allowed */
float lgd_max_cluster_rho;
  /* largest cluster distance from LGD center allowed */
float lgd_min_cluster_separation;
  /* smallest cluster-cluster separation allowed */

int setup_lgd_cluster_cleanup(char *users_filename){
  FILE *config_file=NULL;
  if (users_filename==NULL){
    char *internal_filename;
    char *radphi_config;
    char default_config_filename[MAX_PATH]="/maps/fiducial_config.txt";
    radphi_config=getenv("RADPHI_CONFIG");
    internal_filename=malloc(sizeof(char)*MAX_PATH);
    strcpy(internal_filename,radphi_config);
    strcat(internal_filename,default_config_filename);
#ifdef DEBUGGING
    fprintf(stderr,"users filename: %s\n",internal_filename);
#endif
    if ((config_file=fopen(internal_filename,"r"))==NULL) {
      fprintf(stderr,"Error in setup_lgd_cluster_cleanup:\n");
      fprintf(stderr,"unable to open file %s for configuration.\n",
	      internal_filename);
      free(internal_filename);
      return(1);
    }
    free(internal_filename);
  }
  else {
#ifdef DEBUGGING
    fprintf(stderr,"users filename: %s\n",users_filename);
#endif
    if((config_file=fopen(users_filename,"r"))==NULL){
      fprintf(stderr,"Error in setup_lgd_cluster_cleanup:\n");
      fprintf(stderr,"unable to open file %s for configuration.\n",
	      users_filename);
      return(1);
    }
  }
  fscanf(config_file,"%f",&lgd_min_cluster_rho);
  if(seek_past_next_newline(config_file)){
    fprintf(stderr,"Error in setup_lgd_cluster_cleanup:\n");
    fprintf(stderr,"Error finding next newline in config file.\n");
    fclose(config_file);
    return(1);
  }
#ifdef DEBUGGING
  fprintf(stderr,"got inner fid cut: %.1f cm.\n",lgd_min_cluster_rho);
#endif

  fscanf(config_file,"%f",&lgd_max_cluster_rho);
  if(seek_past_next_newline(config_file)){
    fprintf(stderr,"Error in setup_lgd_cluster_cleanup:\n");
    fprintf(stderr,"Error finding next newline in config file.\n");
    fclose(config_file);
    return(1);
  }
#ifdef DEBUGGING
  fprintf(stderr,"got outer fid cut: %.1f cm.\n",lgd_max_cluster_rho);
#endif

  fscanf(config_file,"%f",&lgd_min_cluster_separation);
  if(feof(config_file)){
    fprintf(stderr,"Error in setup_lgd_cluster_cleanup:\n");
    fprintf(stderr,"unexpected end of config file.\n");
    fclose(config_file);
    return(1);
  }  
#ifdef DEBUGGING
  fprintf(stderr,"got separation cut: %.2f cm.\n",
	  lgd_min_cluster_separation);
#endif
  lgd_cluster_cleanup_initialized=1;
  fclose(config_file);
  return(0);
}

int lgd_cluster_cleanup(itape_header_t *event){
  lgd_clusters_t *clusters=NULL;
  int i,j;
  float rho;
  float sep;
  float comb_width;
  int cleanup_result=0;

  if((clusters=data_getGroup(event,GROUP_LGD_CLUSTERS,0))==NULL){
#ifdef VERBOSE
    fprintf(stderr,"Warning in lgd_cluster_cleanup:\n");
    fprintf(stderr,"called with event with no clusters!\n");
#endif
    return(1);
  }
  for(i=0;i<clusters->nClusters;i++){
    /* i goes over all the cluster indices */
    rho=sqrt((double)(SQR(clusters->cluster[i].space.x)+
	     SQR(clusters->cluster[i].space.y)));
    if(rho<lgd_min_cluster_rho)
      cleanup_result|=CLUSTER_CLEANUP_INNER_EDGE; 
    if(rho>lgd_max_cluster_rho)
      cleanup_result|=CLUSTER_CLEANUP_OUTER_EDGE;
#if defined CLUSTER_CLEANUP_SEPARATION_CUT   
    for(j=i+1;j<clusters->nClusters;j++){
      sep=sqrt(SQR(clusters->cluster[i].space.x-
		   clusters->cluster[j].space.x)+
	       SQR(clusters->cluster[i].space.y-
		   clusters->cluster[j].space.y));
      comb_width=(clusters->cluster[i].width + clusters->cluster[j].width);
/*      if((sep - comb_width) < 10.0) 
	cleanup_result|=CLUSTER_CLEANUP_SEPARATION; */
      if (CLUSTER_CLEANUP_SEPARATION_CUT) {
        if(sep<CLUSTER_CLEANUP_SEPARATION_CUT)
	  cleanup_result|=CLUSTER_CLEANUP_SEPARATION;
      }
      else {
        if(sep<lgd_min_cluster_separation)
	  cleanup_result|=CLUSTER_CLEANUP_SEPARATION;
      }	      
	      
    }
#endif	
  }
  return(cleanup_result);
}
      
int seek_past_next_newline(FILE *filept){
  char byte;
  while((byte=fgetc(filept))!='\n'){
    if(byte==EOF)
      return(1);
  }
  return(0);
}

/* Craig Steffen 19MAY98
   removing one block clusters with < .25 of energy */
    

void prune_cluster_singles(lgd_clusters_t *clusters, lgd_hits_t *clusterHits)
{
   int nhits=0;
   int nClusters=0;
   int ci,hi; /* cluster index, hit index */

   nhits=0;
   nClusters=0;
   for(ci=0;ci<clusters->nClusters;ci++){
      if( (clusters->cluster[ci].nBlocks == 1)
#if defined  SINGLE_BLOCK_CLUSTER_THRESHOLD		    
	  && (clusters->cluster[ci].energy < SINGLE_BLOCK_CLUSTER_THRESHOLD)
#endif
	 ){
	/* rejected cluster stuff here */
      }
      else{
	clusters->cluster[nClusters].nBlocks = clusters->cluster[ci].nBlocks;
	clusters->cluster[nClusters].flags = clusters->cluster[ci].flags;
	clusters->cluster[nClusters].space = clusters->cluster[ci].space;
	clusters->cluster[nClusters].width = clusters->cluster[ci].width;
	clusters->cluster[nClusters].energy = clusters->cluster[ci].energy;
	clusters->cluster[nClusters].firstBlock = nhits;
	for(hi=clusters->cluster[ci].firstBlock;
            hi<clusters->cluster[ci].firstBlock+clusters->cluster[ci].nBlocks;
            hi++) {
	   clusterHits->hit[nhits].channel = clusterHits->hit[hi].channel;
	   clusterHits->hit[nhits].energy = clusterHits->hit[hi].energy;
	   nhits++;
	}
	nClusters++;
      }
   }
   clusters->nClusters = nClusters;
   clusterHits->nhits = nhits;
}

void splitoff_reclaim(itape_header_t *event, lgd_clusters_t *clusters)
{
  int i, j, k, l;
  int split_mult=0;
  int pair[2];
  float sep;
  float lowest_sep=120.0;
  float sum_energy, sum_width;
  float weight_i, weight_j, sum_weight;
  float power;
  float size;
  
  lgd_reconstruct_t *split_clusters, *temp;

  /**** Give intial values to split_clusters ****/
  split_clusters = malloc(sizeof_lgd_reconstruct_t(clusters->nClusters));

  split_clusters->nClusters = clusters->nClusters;
  split_clusters->split_mult = 0;
  split_clusters->next_conf_level = 1.0;
  for(k=0;k<clusters->nClusters;k++) {
    split_clusters->cluster[k].energy = clusters->cluster[k].energy;
    split_clusters->cluster[k].width = clusters->cluster[k].width;
    split_clusters->cluster[k].space.x = clusters->cluster[k].space.x;
    split_clusters->cluster[k].space.y = clusters->cluster[k].space.y;
    split_clusters->cluster[k].space.z = clusters->cluster[k].space.z;
    split_clusters->cluster[k].flags = clusters->cluster[k].flags;
    split_clusters->cluster[k].firstBlock = clusters->cluster[k].firstBlock;
    split_clusters->cluster[k].nBlocks = clusters->cluster[k].nBlocks;
  }

  /**** Iteratively produce reconstructions ****/
  for(split_mult=0;split_mult<clusters->nClusters;split_mult++) {

    if(split_mult != 0) {

      /* Sum clusters */
      sum_energy = split_clusters->cluster[pair[0]].energy + split_clusters->cluster[pair[1]].energy;
      weight_i = log(split_clusters->cluster[pair[0]].energy / sum_energy);
      weight_j = log(split_clusters->cluster[pair[1]].energy / sum_energy);
      sum_weight = weight_i + weight_j;

      split_clusters->cluster[pair[0]].space.x = 
	((weight_i*split_clusters->cluster[pair[0]].space.x)+
	 (weight_j*split_clusters->cluster[pair[1]].space.x)) / sum_weight;
      split_clusters->cluster[pair[0]].space.y = 
	((weight_i*split_clusters->cluster[pair[0]].space.y)+
	 (weight_j*split_clusters->cluster[pair[1]].space.y)) / sum_weight;
      split_clusters->cluster[pair[0]].space.z = 
	((weight_i*split_clusters->cluster[pair[0]].space.z)+
	 (weight_j*split_clusters->cluster[pair[1]].space.z)) / sum_weight;
      split_clusters->cluster[pair[0]].energy = sum_energy;
      split_clusters->cluster[pair[1]].energy = -1.0;

      /* Compress clusters */
      for(l=pair[1];l<(clusters->nClusters-1);l++) {
	split_clusters->cluster[l].energy = split_clusters->cluster[l+1].energy;
	split_clusters->cluster[l].width = split_clusters->cluster[l+1].width;
	split_clusters->cluster[l].space.x = split_clusters->cluster[l+1].space.x;
	split_clusters->cluster[l].space.y = split_clusters->cluster[l+1].space.y;
	split_clusters->cluster[l].space.z = split_clusters->cluster[l+1].space.z;
	split_clusters->cluster[l].flags = split_clusters->cluster[l+1].flags;
	split_clusters->cluster[l].firstBlock = split_clusters->cluster[l+1].firstBlock;
	split_clusters->cluster[l].nBlocks = split_clusters->cluster[l+1].nBlocks;
      }
      split_clusters->nClusters--;
    }

    /** Find next most likely splitoff pair **/
    lowest_sep=120.0;
    for(i=0;i<split_clusters->nClusters;i++) {
      for(j=i+1;j<split_clusters->nClusters;j++) {
	sep=sqrt(SQR(split_clusters->cluster[i].space.x - split_clusters->cluster[j].space.x)+
		 SQR(split_clusters->cluster[i].space.y - split_clusters->cluster[j].space.y));
	if(sep < lowest_sep) {
	  lowest_sep=sep;
	  pair[0]=i; pair[1]=j;
	}
      }
    }

    split_clusters->next_conf_level = lowest_sep;
    split_clusters->split_mult = split_mult;
      
    /* Write to subgroup of GROUP_LGD_RECONSTRUCT */
    size = sizeof_lgd_reconstruct_t(split_clusters->nClusters);
    temp = data_addGroup(event,BUFSIZE,GROUP_LGD_RECONSTRUCT,split_clusters->nClusters,size);
    memcpy(temp,split_clusters,size);
  }
  free(split_clusters);
}

/*** 
     DS 3-18-03
     <apply_reclamation> is a quick measure which goes through the subgroups 
     filled above until there are no more clusters closer than the value stated
     as "min_cluster_sep" - It currently fills events into the GROUP_LGD_CLUSTERS
     to facilitate switching on and off of this tool (without rewriting my 
     version of RODD), but these could be placed somewhere else by someone less 
     confident in the reconstruction
***/

void apply_reclamation(itape_header_t *event, lgd_clusters_t *clusters)
{
  int i,j;
  int nClusters;
  float min_cluster_sep = CLUSTER_RECLAIM_THRESHOLD;
  lgd_reconstruct_t *reconstruct;

  nClusters = clusters->nClusters;
  for(i=0;i<clusters->nClusters;i++){
    reconstruct = data_getGroup(event,GROUP_LGD_RECONSTRUCT,(nClusters - i));
    for(j=0;j<reconstruct->nClusters;j++){
      clusters->cluster[j].energy = reconstruct->cluster[j].energy;
      clusters->cluster[j].width = reconstruct->cluster[j].width;
      clusters->cluster[j].space.x = reconstruct->cluster[j].space.x;
      clusters->cluster[j].space.y = reconstruct->cluster[j].space.y;
      clusters->cluster[j].space.z = reconstruct->cluster[j].space.z;
      clusters->cluster[j].flags = reconstruct->cluster[j].flags;
      clusters->cluster[j].firstBlock = reconstruct->cluster[j].firstBlock;
      clusters->cluster[j].nBlocks = reconstruct->cluster[j].nBlocks;
    }
    clusters->nClusters = reconstruct->nClusters;
    if(reconstruct->next_conf_level > min_cluster_sep) break;
  }
}


int make_lgd_clusters(itape_header_t *event, lgd_clusters_t *p_clust, int max)
{
  lgd_hits_t *lgdHits=NULL;
  lgd_hits_t *clusterHits=NULL;
  lgd_clusters_t *clusters=NULL;
  itape_header_t *tmp;

  int size, ret;

  if (p_clust == NULL) {
    clusters = malloc(sizeof_lgd_hits_t(max));
  }
  else{
    clusters = p_clust;
  }

  lgdHits = data_getGroup(event,GROUP_LGD_HITS,0);

/*   associate hits into clusters  */
  if (lgdHits) {
#if defined LGD_CLUSTERIZER_UC
    lgdClusterizeUC(lgdHits,&clusterHits,&clusters);
#elif defined LGD_CLUSTERIZER_IU
    lgdClusterizeIU(lgdHits,&clusterHits,&clusters); /* buffer overflow safe? */
#else
#   error No clusterizer specified - must select either \
          LGD_CLUSTERIZER_IU or LGD_CLUSTERIZER_UC !
#endif
  }
  else {
    static int first_time=1;
    if (first_time) {
#ifdef VERBOSE
      fprintf(stderr,"Warning from make_lgd_clusters: ");
      fprintf(stderr,"cannot clusterize without lgd hits.\n");
#endif
      first_time=0;
    }
    clusters->nClusters=0;
  }

  /*   prune undesireable clusters  - Craig Steffen */
#if defined PRUNE_SINGLES
  if (clusterHits) {
    prune_cluster_singles(clusters,clusterHits);
  }
#endif 

  /** Reclaim splitoff clusters: DJS
      Fills sub-GROUP_LGD_RECONSTRUCT's, then applies cut and
      writes results to a lgd_cluster_t GROUP **/
  
#if defined SPLITOFF_RECLAMATION  
  splitoff_reclaim(event, clusters);
  apply_reclamation(event, clusters);
#endif  

  if (p_clust == NULL) {
    int size = sizeof_lgd_clusters_t(clusters->nClusters);
    photons_t *tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_CLUSTERS,0,size);
    memcpy(tmp,clusters,size);
    free(clusters); 
  }

  if (lgdHits) {
    int size = sizeof_lgd_hits_t(clusterHits->nhits);
    lgd_clusters_t *tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_CLUSTER_HITS,0,size);
    memcpy(tmp,clusterHits,size);
    free(clusterHits);
  }
  return(0);
}
