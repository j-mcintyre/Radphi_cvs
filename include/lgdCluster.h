/* lgdCluster.h 4.2 
   06 Jun 1996 22:19:43
   /auto/home/lemond/e852/source/include/SCCS/s.lgdCluster.h
*/

#include <ntypes.h>

#ifndef LGDCOORSH_INCLUDED
#define LGDCOORSH_INCLUDED

#define INDIANA_ALGORITHM   1
#define NOTREDAME_ALGORITHM 2
#define LGDCLUS_LINEAR 1
#define LGDCLUS_LOG 2
#define LGDCLUS_LOG2 3
#define LGDCLUS_BUGGE 4 /* Not recommended */

#define BUGGE_MAGIC_THRESHOLD 0.01
#define BUGGE_DEFAULT_PARAM 0.975

/* LGD clusterizer default parameters changed
   to be same as listed in 'http://anthrax.physics.indiana.edu/Misc-E852/standardCuts.html'
   C.O. 1996-May-12
*/

#define LGDCLUS_DEF_LEVEL          3
#define LGDCLUS_DEF_SEEDENERGY1    0.35
#define LGDCLUS_DEF_SEEDENERGY2    0.15
#define LGDCLUS_DEF_POSITION       LGDCLUS_LOG2
#define LGDCLUS_DEF_LOGOFFSET      4.1  /* not used by LGDCLUS_LOG2 */
#define LGDCLUS_DEF_MAXNUMCLUSTERS 8
#define LGDCLUS_DEF_MAXNUMSHARED   3
#define LGDCLUS_DEF_MAXRASSOCIATE 10.0
#define LGDCLUS_DEF_SINGLEBLOCKWIDTH 1.155
#define LGDCLUS_DEF_MAXRCHARGED   10.0

#define LGDCLUS_TRUE 1
#define LGDCLUS_FALSE 0

typedef struct {
  int channel;
  float energy;
} lgd_hit_t;


typedef struct
{
  int nhits;
  lgd_hit_t hit[1];
}lgd_hits_t;

#define sizeof_lgd_hits_t(N) \
       (sizeof(lgd_hits_t) + (((N)>0)?(N)-1:0)*sizeof(lgd_hit_t))

#define LGD_F_HADRONIC    0x0010
#define LGD_F_DESTREAKED  0x0020

typedef struct
{
  int nBlocks;          /* number of blocks in this cluster            */
  int firstBlock;       /* index into the GROUP_LGD_CLUSTER_HITS group */
  int flags;            /* LGD_F_xxx flags defined above               */
  vector3_t space;      /* cluster position in the LGD local position  */
  float width;          /* width of the cluster in 'cm'                */
  float energy;         /* energy of the cluster in 'GeV'              */
}lgd_cluster_t;

typedef struct{
  int nClusters;
  lgd_cluster_t cluster[1];
}lgd_clusters_t;

#define sizeof_lgd_clusters_t(N) \
       (sizeof(lgd_clusters_t) + (((N)>0)?(N)-1:0)*sizeof(lgd_cluster_t))

typedef struct{
  int nClusters;
  int split_mult;
  float next_conf_level;
  lgd_cluster_t cluster[1];
}lgd_reconstruct_t;

#define sizeof_lgd_reconstruct_t(N) \
       (2*sizeof(int) + sizeof(float) + (N)*sizeof(lgd_cluster_t))

void lgdSetupAdvanced(int,float,float,int,float,int,int,float,float);
void lgdSetup(int,float,float,int,float);
int lgdClusterizeIU(lgd_hits_t *,lgd_hits_t **,lgd_clusters_t **);
int lgdClusterizeUC(lgd_hits_t *,lgd_hits_t **,lgd_clusters_t **);

#endif
