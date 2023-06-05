/*
 * iitypes.h - unpacked data after gain and pedestal corrections
 *
 * The unpacking is logically divided into stage i and stage ii, with
 *
 *    stage i  : raw adc and tdc value structures independent of the
 *               meaning of the value or the units of the quantity
 *               it measures, defined in itypes.h
 *
 *    stage ii : pulse-height and time information structures on a
 *               per-detector basis expressed as floats in physical
 *               units of GeV (deposited energy) or ns (time since
 *               t0 for the event), defined in iitypes.h
 *
 * The stage ii structures for the lgd are defined in lgdCluster.h
 * but follow the same spirit as those defined below for the RPP,
 * CPV, UPV and the tagger.
*/

#ifndef IITYPES_INCLUDED
#define IITYPES_INCLUDED


#include <ntypes.h>          /* for int32, etc... definitions */
#include <particleType.h>    /* for particle type codes       */

/*
 * The time_list structure is set up as a shared structure for the tdc
 * hit information from all of the detectors.  In order to be able to
 * read/write the structures between a disk file and memory, one should
 * avoid using pointers and use array indices instead.  The time_list is
 * a single array of time values for all hits, with hits for a single
 * counter appearing contiguously.  The individual detector hit structure
 * contains a count of tdc hits for that counter followed by the index
 * into time_list of the first hit in the sequence.  For generality, even
 * the detectors being read out with single-hit tdcs use the time_list
 * structures with one element per counter.  Set up this way, it is easy
 * to know the size needed for the time_list structure before building it.
*/

typedef struct
{
  float32 le;		/* time since t0 (ns) */
} time_list_el_t;

typedef struct
{
  uint32 nelements;
  time_list_el_t element[1];
} time_list_t;

#define sizeof_time_list_t(N) \
       (sizeof(time_list_t) + (((N)>0)?(N)-1:0)*sizeof(time_list_el_t))

typedef struct
{
  int32 channel;
  float32 energy;	/* GeV */
  uint32 times;
  uint32 tindex;
} cpv_hit_t;

typedef struct
{
  uint32 nhits;
  cpv_hit_t hit[1];
} cpv_hits_t;

#define sizeof_cpv_hits_t(N) \
       (sizeof(cpv_hits_t) + (((N)>0)?(N)-1:0)*sizeof(cpv_hit_t))

typedef struct
{
  int32 channel;
  float32 energy;	/* GeV */
  uint32 times;
  uint32 tindex;
} upv_hit_t;

typedef struct
{
  uint32 nhits;
  upv_hit_t hit[1];
} upv_hits_t;

#define sizeof_upv_hits_t(N) \
       (sizeof(upv_hits_t) + (((N)>0)?(N)-1:0)*sizeof(upv_hit_t))

typedef struct
{
  int32 channel;
  float32 energy;
  uint32 times;
  uint32 tindex;
} bsd_hit_t;

typedef struct
{
  uint32 nhits;
  bsd_hit_t hit[1];
} bsd_hits_t;

#define sizeof_bsd_hits_t(N) \
       (sizeof(bsd_hits_t) + (((N)>0)?(N)-1:0)*sizeof(bsd_hit_t))

typedef struct
{
  int32 channel;
  float32 energy[2];    /*  order is upstream, downstream  */
  uint16 times[2];
  uint16 tindex[2];
} bgv_hit_t;

typedef struct
{
  uint32 nhits;
  bgv_hit_t hit[1];
} bgv_hits_t;

#define sizeof_bgv_hits_t(N) \
       (sizeof(bgv_hits_t) + (((N)>0)?(N)-1:0)*sizeof(bgv_hit_t))

typedef struct
{
  int32 channel;
  float32 energy;	/* GeV */
  uint32 times;
  uint32 tindex;
} rpd_hit_t;

typedef struct
{
  uint32 nhits;
  rpd_hit_t hit[1];
} rpd_hits_t;

#define sizeof_rpd_hits_t(N) \
       (sizeof(rpd_hits_t) + (((N)>0)?(N)-1:0)*sizeof(rpd_hit_t))

/*
 * The following structures supplement the above RPD data
 * to provide ready information on G/E counter coincidences
*/

typedef struct
{
  uint32 Ehit;
  uint32 Ghit;
} rpd_coin_t;

typedef struct
{
  uint32 ncoins;
  rpd_coin_t coin[1];
} rpd_coins_t;

#define sizeof_rpd_coins_t(N) \
       (sizeof(rpd_coins_t) + (((N)>0)?(N)-1:0)*sizeof(rpd_coin_t))

/*
 * The following structures supplement the above BSD data
 * to provide ready information on R,L,S pixels
*/

typedef struct
{
  uint32 right;
  uint32 left;
  uint32 straight;
  uint32 ring;
  float  z;
  float  phi;
  float energy[3]; /* [0] is energy_r, [1] is energy_l, [2] is energy_s*/
  float time[3]; /* [0] is time_r, [1] is time_l, [2] is time_s */
} bsd_pixel_t;

typedef struct
{
  uint32 npixels;
  bsd_pixel_t pixel[1];
} bsd_pixels_t;

#define sizeof_bsd_pixels_t(N) \
       (sizeof(bsd_pixels_t) + (((N)>0)?(N)-1:0)*sizeof(bsd_pixel_t))


/* The following is the structure type that will be used to hold the 
 * recoil information that can be found from the barrel detectors.
 * - npixels      : number of BSD pixels that make up this recoil cluster
 * - pixel_index  : an array that contains indices (corresponding to the
 *                  pixels group) of the pixels that made up the cluster,
 *                  so that they can be found to check time, add up energy,
 *                  and that sort of thing.  
 * - theta,phi    : polar angle (from the +Z axis) and azimuthal angle
 *                  of the cluster centroid (radians)
 * - time         : average time-of-impact for recoil particle on the
 *                  barrel scintillators, with light propagation delay
 *                  removed (ns).
 * - dEdx         : average dE/dx of particle in barrel scintillators,
 *                  after path-length correction (GeV/cm).
 * - energy       : total energy of recoil particle, based on the summed
 *                  energy deposition in the BSD and BGD counters (GeV).
 */ 

#define MAX_PIXELS_IN_RECOIL 5

typedef struct{
  uint32 npixels;
  uint32 pixel_index[MAX_PIXELS_IN_RECOIL];
  float theta;
  float phi;
  float time;
  float dEdx;
  float energy;
} recoil_t;

typedef struct{
  uint32 nrecoils;
  float32 best_t0;	/* best estimate for event t0, based on pixel times */
  recoil_t recoil[1];
} recoils_t;

#define sizeof_recoils_t(N) \
       (sizeof(recoils_t) + (((N)>0)?(N)-1:0)*sizeof(recoil_t))

typedef struct
{
  uint32 channel;
  float32 energy[2];	/* limits [lower,upper] of tagging bin (GeV) */
  uint32 times[2];	/* first left, then right */
  uint32 tindex[2];	/* first left, then right */
} tagger_hit_t;

typedef struct
{
  uint32 nhits;
  tagger_hit_t hit[1];
} tagger_hits_t;

#define sizeof_tagger_hits_t(N) \
       (sizeof(tagger_hits_t) + (((N)>0)?(N)-1:0)*sizeof(tagger_hit_t))

typedef struct
{
  uint32 channel;
  float32 energy;
  float32 time;
  uint32 next;
} tagger_cluster_t;

typedef struct
{
  int nclusters;
  tagger_cluster_t cluster[1];
} tagger_clusters_t;
  
#define sizeof_tagger_clusters_t(N) \
       (sizeof(tagger_clusters_t) + (((N)>0)?(N)-1:0)*sizeof(tagger_cluster_t))

typedef struct
{
  uint32 channel;
  float32 energy;
  float32 time;
  float32 tagging_weight;
} tagger_photon_t;

typedef struct
{
  int nphotons;
  tagger_photon_t photon[1];
} tagger_photons_t;
  
#define sizeof_tagger_photons_t(N) \
       (sizeof(tagger_photons_t) + (((N)>0)?(N)-1:0)*sizeof(tagger_photon_t))

typedef struct
{
  int cluster;          /* index in cluster table, <0 if barrel cluster   */
  vector3_t momentum;   /* px, py, pz in 'GeV/c'                          */
  vector3_t position;   /* cluster position at detector in 'cm '          */
  float energy;         /* energy of the photon in 'GeV'                  */
} photon_t;

typedef struct
{
  int nPhotons;          /* number of photons   */
  photon_t photon[1];   
} photons_t;

#define sizeof_photons_t(N) \
       (sizeof(photons_t) + (((N)>0)?(N)-1:0)*sizeof(photon_t))

typedef struct
{
  Particle_t type;         /* type of meson, defined in particleType.h    */
  vector3_t momentum;      /* px, py, pz in 'GeV/c'                       */
  vector3_t vertex;        /* decay vertex position at detector in 'cm '  */
  float energy;            /* energy of the meson in 'GeV'                */
  float mass;              /* mass of meson in 'GeV'                      */
  int ndecay;              /* number of decay products, if any, <= 3      */
  Particle_t dtype[3];     /* particle type of daughter particles         */
  int daughter[3];         /* index of meson or photon daughters into     */
} meson_t;                 /*         mesons_t table or photons_t table   */

typedef struct
{
  int nMesons;             /* number of mesons */
  meson_t meson[1];
} mesons_t;

#define sizeof_mesons_t(N) \
       (sizeof(mesons_t) + (((N)>0)?(N)-1:0)*sizeof(meson_t))


#define MAX_HYPOTHESIS_FITS 100

#define FIT_TYPE_UNINITIALIZED                -1
#define FIT_TYPE_DUMMY_FIT                     1

#define FIT_TYPE_THREE_PHOTON_NO_HYPOTHESIS  300

#define FIT_TYPE_THREE_PHOTON_PI0_GAMMA_12_0 310
#define FIT_TYPE_THREE_PHOTON_PI0_GAMMA_02_1 311
#define FIT_TYPE_THREE_PHOTON_PI0_GAMMA_01_2 312

#define FIT_TYPE_THREE_PHOTON_ETA_GAMMA_12_0 320
#define FIT_TYPE_THREE_PHOTON_ETA_GAMMA_02_1 321
#define FIT_TYPE_THREE_PHOTON_ETA_GAMMA_01_2 322

/* photon fit structures */
typedef struct{
  int type;
  vector3_t photon[3];
  float chi_sq;
  float penalty_func;
  float constraint_func;
  float beam_E;
  recoils_t recoil;
} three_photon_fit_t;

typedef struct{
  int nfits;
  three_photon_fit_t fit[1];
} three_photon_fits_t;

#define sizeof_three_photon_fits_t(N) \
       (sizeof(three_photon_fits_t) + (((N)>0)?(N)-1:0)*sizeof(three_photon_fit_t))

#endif

/* end file */
