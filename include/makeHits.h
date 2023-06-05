#ifndef MAKEHITS_H_INCLUDED
#define MAKEHITS_H_INCLUDED

/* 
 * Instructions for using the makeHits library
 * 1. Before calling any makeHits functions, do a setup_makeHits.  This
 *    needs to be repeated at the start of each new run.  setup_makeHits
 *    does a call to umap_loadDb() itself, so the user does not need to.
 *
 * 2. The simplest way to make hits is to invoke the function makeHits()
 *    that calls the functions for the individual components, and appends
 *    the groups to the end of the existing event in memory.  It assumes
 *    a (large) event buffer of size BUFSIZE.
 *
 * 3. If you call the makers individually then you can put the hits into
 *    buffers separate from the exiting event in memory.  To do this, 
 *    set up a pointer of the appropriate type to point to an area
 *    with plenty of space before passing it to any makeHits function.
 *    If the call is made with the pointer null, the hits group will be
 *    appended to the end of the existing event, as with makeHits.
 *
 * 4. The timelist group is shared by several subsystems and so it is
 *    created automatically when needed, with enough capacity to hold
 *    all of the hits in the event.  It can be accessed by the user 
 *    through the pointer p_timelist below at any time, or written to
 *    the event by calling make_timelist with a null pointer.  This
 *    call comes after all subsystems that can make timelist hits, or
 *    else the timelist recorded in the event will be incomplete.
 *
 * 5. In the "max" argument you give the maximum number of entries that
 *    the buffer can hold.  You can use the sizeof_bsd_hits() (or other
 *    subsystem) to get the max of the buffer for the call to malloc.
 */

#include <ntypes.h>
#include <disData.h>
#include <itypes.h> 
#include <iitypes.h>
#include <eventType.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <libBGV.h>

#define BUFSIZE 100000
#define TMIN_BSD 1.
#define TMIN_BGV 1.
#define TMIN_CPV 1.
#define TMIN_UPV 1.
#define TMIN_TAG 1
#define TMAX_BSD 2200.
#define TMAX_BGV 4100.
#define TMAX_CPV 2200.
#define TMAX_UPV 2200.
#define TMAX_TAG 4100
#define AMIN_BSD 1.
#define AMIN_BGV 1.
#define AMIN_CPV 1.
#define AMIN_UPV 1.
#define AMIN_LGD 1.
#define AMIN_PIXEL 100.
#define TDC_1877_GAIN -0.5

/* default values in LGD nonlinearity correction */
#define BETA 0.010
#define EPSILON 0.020

#define MAX_TIME_LIST_LENGTH    2000
#define MAX_TAGGER_CLUSTERS      300
#define MAX_TAGGER_PHOTONS       300

/* define channel number variables; one assumes that the sizes
   of things will not change as a function of run in this 
   experiment.  CPS June 26, 2000 */

#define N_TAGGER_ENERGY_CHANNELS 19

#define CLUSTER_CLEANUP_INNER_EDGE 0x1
#define CLUSTER_CLEANUP_OUTER_EDGE 0x2
#define CLUSTER_CLEANUP_SEPARATION 0x4

#define N_BSD_PIXELS 192

#define N_CHANNEL_STATUS_ARRAYS 6

#define N_TAG_CHANNELS 38
#define N_TAG_COUNTERS 19
#define N_UPV_CHANNELS 8
#define N_UPV_COUNTERS 8
#define N_BSD_CHANNELS 48
#define N_BSD_COUNTERS 48
#define N_BGV_CHANNELS 48
#define N_BGV_COUNTERS 24
#define N_CPV_CHANNELS 30
#define N_CPV_COUNTERS 30

#define LGD_STATUS_BITMASK_DEFAULT 0x0
#define UNSIGNED_CHAR_ALL_BITS_SET 0xff

#define N_GLOBAL_DETECTOR_NUMBERS 6
#define GLOBAL_DETECTOR_NUMBER_TAG 0
#define GLOBAL_DETECTOR_NUMBER_UPV 1
#define GLOBAL_DETECTOR_NUMBER_BSD 2
#define GLOBAL_DETECTOR_NUMBER_BGV 3
#define GLOBAL_DETECTOR_NUMBER_CPV 4
#define GLOBAL_DETECTOR_NUMBER_LGD 5

/* hard coded values (these things are not a function of run number,
   but are strictly static for the 2000 data run. */

#define BSD_N_RINGS 8
#define BSD_N_PIXELS_PER_RING 24
#define BSD_N_RIGHT_PADDLES 12
#define BSD_N_LEFT_PADDLES 12
#define BSD_N_STRAIGHT_PADDLES 24
#define TWO_PI (M_PI*2.0)
#define BSD_PHI_ZERO (.392699081699)
#define BSD_PHI_INCREMENT (.261799387799)
#define BSD_RADIUS 35.0
#define TAG_BSD_TIME_WINDOW_HALF_WIDTH (5.0) /* (this is in ns) */

#define TAG_LR_MAX_DIFFERENCE 1.75 /* this is in ns */

typedef struct {
   Particle_t meson;
   float width;
   Particle_t daughter[2];
} meson_search_list_el_t;

typedef struct {
   int stages;
   meson_search_list_el_t try[1];
} meson_search_list_t;

#define sizeof_meson_search_list_t(N) \
       (sizeof(meson_search_list_t) + (((N)>0)?(N)-1:0)*sizeof(meson_search_list_el_t))

extern int monte_carlo;
extern int thisRun;
extern int thisEvent;
extern float *randomNo;
extern int randomNos, randomNext;
extern float *bgv_ped, *bgv_cc; extern int *bgv_thresh;
extern float *cpv_ped, *cpv_cc; extern int *cpv_thresh;
extern float *upv_ped, *upv_cc; extern int *upv_thresh;
extern float *bsd_ped, *bsd_cc; extern int *bsd_thresh;
extern float *lgd_ped, *lgd_cc; extern int *lgd_thresh;
extern float *lgd_cc_tune;
extern float *bgv_tdc_offset;
extern float *bgv_offset_ave, *bgv_mean_diff;   // no longer used [rtj]
extern float *bgv_US_offset, *bgv_DS_offset;    // no longer used [rtj]
extern float *bgv_atten_length, *bgv_tdc_offset;
extern float *bsd_atten_length, *bsd_tdc_offset;
extern float bsd_t0_offset;
extern float bgv_t0_offset;
extern float cpv_t0_offset;
extern float upv_t0_offset;
extern float tag_t0_offset;
extern int *bsd_tdc_upper_bound,*bsd_tdc_lower_bound;
extern float *tag_tdc_offset;
extern float *cpv_tdc_offset;
extern float *upv_tdc_offset;
extern int lgd_nChannels, lgd_nRows, lgd_nCols;
extern time_list_t *p_timelist;
extern float epsilon;
extern float beta;
extern float bsd_minus_tagger_offset;
/* channel status definitions */
extern unsigned char *channel_status_tag;
extern unsigned char *channel_status_upv;
extern unsigned char *channel_status_bsd;
extern unsigned char *channel_status_bgv;
extern unsigned char *channel_status_cpv;
extern unsigned char *channel_status_lgd;

extern float bsd_paddle_diff_center_by_pixel_0[N_BSD_PIXELS];
extern float bsd_paddle_diff_width_by_pixel_0[N_BSD_PIXELS];
extern float bsd_paddle_diff_center_by_pixel_1[N_BSD_PIXELS];
extern float bsd_paddle_diff_width_by_pixel_1[N_BSD_PIXELS];
extern float bsd_paddle_diff_center_by_pixel_2[N_BSD_PIXELS];
extern float bsd_paddle_diff_width_by_pixel_2[N_BSD_PIXELS];

/* global utility arrays */
extern int n_channels_in_detector[N_GLOBAL_DETECTOR_NUMBERS];
extern char *detector_name[N_GLOBAL_DETECTOR_NUMBERS];

extern meson_search_list_t *meson_list;

/* tagger_energy_settings[0] = photon beam endpoint energy (GeV) */
/* tagger_energy_settings[1] = lower tagging limit (fraction of endpoint) */
/* tagger_energy_settings[2] = upper tagging limit (fraction of endpoint) */
extern float tagger_energy_settings[3];

/* tagger_fraction_by_channel is the fraction at mid-channel */
extern float tagger_fraction_by_channel[N_TAGGER_ENERGY_CHANNELS];

/* tagger_CPleakage_by_channel is the CP leakage fraction for each channel */
extern float tagger_CPleakage_by_channel[N_TAGGER_ENERGY_CHANNELS];

/* coordinate of center of pixels in global reference frame */
extern float bsd_ring_midz[BSD_N_RINGS];

/* gain correction factors for each pixel */
extern float bsd_pixel_gain[BSD_N_RINGS][N_BSD_COUNTERS];

#ifdef INITIALIZE_HITS_GLOBALS

int monte_carlo=0;
int thisRun=-1;
int thisEvent=-1;
float *randomNo=NULL;
int randomNos=100000, randomNext=0;
float *bgv_ped=NULL, *bgv_cc=NULL; int *bgv_thresh=NULL;
float *cpv_ped=NULL, *cpv_cc=NULL; int *cpv_thresh=NULL;
float *upv_ped=NULL, *upv_cc=NULL; int *upv_thresh=NULL;
float *bsd_ped=NULL, *bsd_cc=NULL; int *bsd_thresh=NULL;
int *bsd_tdc_upper_bound,*bsd_tdc_lower_bound;
float *lgd_ped=NULL, *lgd_cc=NULL; int *lgd_thresh=NULL;
float *lgd_cc_tune=NULL;
float *bgv_offset_ave=NULL, *bgv_mean_diff=NULL;     // no longer used [rtj]
float *bgv_US_offset=NULL, *bgv_DS_offset=NULL;      // no longer user [rtj]
float *bgv_atten_length=NULL, *bgv_tdc_offset=NULL;
float *bsd_atten_length=NULL, *bsd_tdc_offset=NULL;
float *tag_tdc_offset=NULL;
float *cpv_tdc_offset=NULL;
float *upv_tdc_offset=NULL;
float bsd_t0_offset=0;
float cpv_t0_offset=0;
float upv_t0_offset=0;
float bgv_t0_offset=0;
float tag_t0_offset=0;
int lgd_nChannels=-1, lgd_nRows=-1, lgd_nCols=-1;
time_list_t *p_timelist=NULL;
float epsilon=0;
float beta=0;
float bsd_minus_tagger_offset;
unsigned char *channel_status_tag=NULL;
unsigned char *channel_status_upv=NULL;
unsigned char *channel_status_bsd=NULL;
unsigned char *channel_status_bgv=NULL;
unsigned char *channel_status_cpv=NULL;
unsigned char *channel_status_lgd=NULL;

float bsd_paddle_diff_center_by_pixel_0[N_BSD_PIXELS];
float bsd_paddle_diff_width_by_pixel_0[N_BSD_PIXELS];
float bsd_paddle_diff_center_by_pixel_1[N_BSD_PIXELS];
float bsd_paddle_diff_width_by_pixel_1[N_BSD_PIXELS];
float bsd_paddle_diff_center_by_pixel_2[N_BSD_PIXELS];
float bsd_paddle_diff_width_by_pixel_2[N_BSD_PIXELS];

int n_channels_in_detector[N_GLOBAL_DETECTOR_NUMBERS];
char *detector_name[N_GLOBAL_DETECTOR_NUMBERS]=
            {NULL,NULL,NULL,NULL,NULL,NULL};

meson_search_list_t *meson_list=NULL;

float tagger_fraction_by_channel[N_TAGGER_ENERGY_CHANNELS];
float tagger_CPleakage_by_channel[N_TAGGER_ENERGY_CHANNELS];
float tagger_energy_settings[3];
float bsd_ring_midz[BSD_N_RINGS];
float bsd_pixel_gain[BSD_N_RINGS][N_BSD_COUNTERS];

#endif

#define RANDOM randomNo[(randomNext++) % randomNos]-0.5

void setupNonLinearity(float new_beta, float new_epsilon);
int setup_makeHits(itape_header_t *event);
int makeHits(itape_header_t *event);
int make_bsd_hits(itape_header_t *event, bsd_hits_t *p_bsd, int max);
int make_bsd_pixels(itape_header_t *event, bsd_pixels_t *p_pix, int max);
int make_recoil(itape_header_t *event, recoils_t *p_recoils, int maxrecoils);
int make_bgv_hits(itape_header_t *event, bgv_hits_t *p_bgv, int max);
int make_cpv_hits(itape_header_t *event, cpv_hits_t *p_cpv, int max);
int make_upv_hits(itape_header_t *event, upv_hits_t *p_upv, int max);
int make_lgd_hits(itape_header_t *event, lgd_hits_t *p_lgd, int max);
int make_tagger_hits(itape_header_t *event, tagger_hits_t *p_tag, int max);
int make_tagger_clusters(itape_header_t *event, tagger_clusters_t *p_clusters,
			 int max_clusters);
int make_tagger_photons(itape_header_t *event, tagger_photons_t *p_photons,
                        int max_photons);
int make_time_list(itape_header_t *event, time_list_t *p_tlist, int max);
int make_lgd_clusters(itape_header_t *event, lgd_clusters_t *p_clust, int max);
int make_bgv_clusters(itape_header_t *event, bgv_clusters_t *p_clust, int max);
int make_photons(itape_header_t *event, photons_t *p_phot, int max);
int make_mesons(itape_header_t *event, mesons_t *p_mesons, int max);
int setup_makePhotons(itape_header_t *event);
int setup_makeMesons(itape_header_t *event, meson_search_list_t *list);
vector4_t makePhotonPatch(lgd_cluster_t *cluster,vector3_t vertex);
int setup_lgd_cluster_cleanup(char *users_filename);
int lgd_cluster_cleanup(itape_header_t *event);
void set_lgd_status_bitmask(unsigned char bitmask);
unsigned char lgd_status_bitmask(void);
int suppress_dead_channels(itape_header_t *event);

#endif
