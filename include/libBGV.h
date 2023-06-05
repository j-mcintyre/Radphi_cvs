/* Written 7/99 by Tom Bogue
   Finds clusters of hits in the BGV
   A cluster is defined as the largest set of blocks  such that:
    1)Every block in the cluster is adjacent to another in the cluster
    2)There is at least 1 block in the cluster w/ energy > Emin and w/ an upstream 
      and a downstream TDC that give a valid value for z (-10<z<LENGTH_BGV+10)
    3)There are no blocks w/ E<Emin2
    4)The z-position of any block is w/i 2 sigmas of the mean z
	(the mean z is a running mean of blocks already in the cluster)
*/

#ifndef BGVCLUSTERIZE
#define BGVCLUSTERIZE

#include <iitypes.h>
#include <disData.h>
#include <ntypes.h>

typedef struct
{
  float32 energy;	/* corrected energy of shower (GeV) */
  float32 energy_up;	/* uncorrected sum of upstream energy (GeV) */
  float32 energy_dn;	/* uncorrected sum of downstream energy (GeV) */
  float32 phi;		/* azimuth of shower (radians) */
  float32 z;  		/* z coordinate of shower in global coordinates */
  float32 t;		/* mean time of shower (ns) */
  float32 scintEnergy;	/* total energy in adjacent BSD paddles, outer layer */
} bgv_cluster_t;

typedef struct
{
  int32 nClusters;
  bgv_cluster_t clusters[1];
} bgv_clusters_t;

#define sizeof_bgv_clusters_t(N) \
       (sizeof(bgv_clusters_t) + (((N)>0)?(N)-1:0)*sizeof(bgv_cluster_t))

void bgvSetup(float Emin,float Emin2); /* sets up the system */

/* finds the clusters in hits and puts them in a newly created structure */
/* void bgvClusterize(itape_header_t *event, bgv_clusters_t **clusters);*/
void bgvClusterize(itape_header_t *event,bgv_clusters_t **clusters);

/* returns the 4-vector of the photon that created cluster */
vector4_t makePhotonBGV(const bgv_cluster_t *cluster);

/* smears the MC from the BGV by experimental resolution */
void MCsmearBGV(bgv_hits_t *hits, time_list_t *tlist);

#endif
