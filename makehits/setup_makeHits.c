/*     setup_makeHits.c     */
/*     D.S. Armstrong, R.T. Jones    June 9 2000  */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include <itypes.h> 
#include <eventType.h>
#include <ntypes.h> 

#include <umap.h>
#include <disData.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <lgdGeom.h>
#include <libBGV.h>
#include <detector_params.h>
#include <pedestals.h>
#include <calibration.h> 
#include <detector_params.h>

#define INITIALIZE_HITS_GLOBALS 1

#include <makeHits.h>
/* #define DEBUGGING 1 */


/*   this code sets up the LGD pedestals and thresholds and 
     the clusterizing parameters. It needs to be called at the
     start of the analysis, and each time the run number changes.
     
     a call to it could look like:

        if (setup_makeHits((itape_header_t *) event)){
	  fprintf(stderr,"Problem with setup_makeHits; exiting");
	  fprintf(logfi,"Problem with setup_makeHits; exiting");
	  return;
         }      
*/


int setup_makeHits(itape_header_t *event)
{
  int actualRun;
  int actualRunList[N_CHANNEL_STATUS_ARRAYS];
  int i,j,ch;

  if (lgdGeomSetup(event->runNo))
    return(1);
  lgdGetGeomParam(&lgd_nChannels,&lgd_nRows,&lgd_nCols);

  if (setup_global_utility_structures()){
    fprintf(stderr,"Error setting up global utility structures\n");
    return(1);
  }

  if (data_getGroup(event,GROUP_MC_EVENT,0) !=NULL){
    printf("Analysing Monte Carlo events...\n");
    monte_carlo = 1;
  }
  
  if (umap_loadDb(event->runNo)) {
    fprintf(stderr,"Error in setup_makeHits: database access failure\n");
    return(-1);
  }

  if (umap_loadDb(event->runNo)) {
    fprintf(stderr,"Second database access failed\n");
    return(-1);
  }

  thisRun = event->runNo;
  thisEvent = -1;

  /* Load pedestals, calibration constants, thresholds for BSD  */

  if (bsd_ped)
    free(bsd_ped);
  bsd_ped = malloc(N_BSD_CHANNELS*sizeof(float));
  if (bsd_thresh)
    free(bsd_thresh);
  bsd_thresh = malloc(N_BSD_CHANNELS*sizeof(int));
  if (pedestalLoad(PED_BSD,event->runNo,
                   N_BSD_CHANNELS,bsd_ped,bsd_thresh,&actualRun))
    return(1);

  if (bsd_cc)
    free(bsd_cc);
  bsd_cc = malloc(N_BSD_CHANNELS*sizeof(float));
  if (ccLoad(CC_BSD,event->runNo,N_BSD_CHANNELS,bsd_cc,&actualRun))
    return(1);

  if (bsd_tdc_offset)
    free(bsd_tdc_offset);
  bsd_tdc_offset = malloc(N_BSD_CHANNELS*sizeof(float));
  if (tdcOffsetLoad(TDC_OFFSET_BSD,event->runNo,N_BSD_CHANNELS,
		     bsd_tdc_offset,&bsd_t0_offset,&actualRun)){
    fprintf(stderr,
	      "Error reading BSD TDC offset run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "tdcOffsetLoad (BSD) success: run %d, real run=%d\n",
	      event->runNo,actualRun);
  }
#endif

  if (bsd_tdc_upper_bound)
    free(bsd_tdc_upper_bound);
  bsd_tdc_upper_bound=malloc(N_BSD_CHANNELS*sizeof(float));
  if (tdcBoundLoad(TDC_UPPER_BOUND_BSD,event->runNo,N_BSD_CHANNELS,
			  bsd_tdc_upper_bound,&actualRun)){
    fprintf(stderr,
	      "Error: get BSD tdc upper bounds run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "tdcBoundLoad for BSD upper success: run %d, real run=%d.\n",
	      event->runNo,actualRun);
  }
#endif

  if (bsd_tdc_lower_bound)
    free(bsd_tdc_lower_bound);
  bsd_tdc_lower_bound=malloc(N_BSD_CHANNELS*sizeof(float));
  if (tdcBoundLoad(TDC_LOWER_BOUND_BSD,event->runNo,N_BSD_CHANNELS,
		     bsd_tdc_lower_bound,&actualRun)){
    fprintf(stderr,
	      "Error: get BSD tdc lower bound run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "tdcBoundLoad for BSD lower success: run %d, real run=%d.\n",
	      event->runNo,actualRun);
  }
#endif

#if VERBOSE
  fprintf(stderr,"Loading BSD tagger time offset and pixel difference times.\n");
#endif
  if (tdcOffsetLoad(TIME_OFFSET_BSD_TAG,event->runNo,1,
	 	     &bsd_minus_tagger_offset,&bsd_t0_offset,&actualRun)) {
    fprintf(stderr,
	      "Error reading BSD to TAG offset run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "tdcOffsetLoad (BSD-TAG offset) load success: run %d, real run=%d, offset= %+.2f\n",
	      event->runNo,actualRun,bsd_minus_tagger_offset);
  }
#endif

#if VERBOSE
  fprintf(stderr,"Loading BGV attenuation constants.\n");
#endif
  if (bgv_atten_length)
    free(bgv_atten_length);
  bgv_atten_length = malloc(N_BGV_CHANNELS*sizeof(float));
  if (attenLengthLoadBGV(event->runNo,N_BGV_CHANNELS,
			  bgv_atten_length,&actualRun)){
    fprintf(stderr,
	      "Error: get BGV atten length run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "attenLengthLoadBGV success: run %d, real run=%d.\n",
	      event->runNo,actualRun);
  }
#endif

#if VERBOSE
  fprintf(stderr,"Loading BSD attenuation constants.\n");
#endif
  if (bsd_atten_length)
    free(bsd_atten_length);
  bsd_atten_length = malloc(N_BSD_CHANNELS*sizeof(float));
  if (attenLengthLoadBSD(event->runNo,N_BSD_CHANNELS,
			  bsd_atten_length,&actualRun)){
    fprintf(stderr,
	      "Error: get BSD atten length run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "attenLengthLoadBSD success: run %d, real run=%d.\n",
	      event->runNo,actualRun);
  }
#endif

  /* Load pedestals, calibration constants, thresholds for LGD  */

  if (lgd_ped)
    free(lgd_ped);
  lgd_ped = malloc(lgd_nChannels *sizeof(float));
  if (lgd_thresh)
    free(lgd_thresh);
  lgd_thresh = malloc(lgd_nChannels *sizeof(int));
  if (pedestalLoad(PED_LGD,event->runNo,lgd_nChannels,
		  lgd_ped,lgd_thresh,&actualRun)){
    fprintf(stderr,"pedestalLoad FAILED!  run=%d\n",event->runNo);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	  "pedestalLoad successful for run %d, with %d channels, real run=%d.\n",
	  event->runNo,lgd_nChannels,actualRun);
  }
#endif

  if (lgd_cc)
    free(lgd_cc);
  lgd_cc = malloc(lgd_nChannels *sizeof(float));
  if (ccLoad(CC_LGD,event->runNo,lgd_nChannels,lgd_cc,&actualRun)){
    fprintf(stderr,"Error reading calibration constants for run %d!\n",event->runNo);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	  "ccLoad successful for run %d, real run=%d.\n",
	  event->runNo,actualRun);
  }
#endif

  if (lgd_cc_tune)
    free(lgd_cc_tune);
  lgd_cc_tune = malloc(lgd_nChannels *sizeof(float));
  if (lgdTuneLoad(event->runNo,lgd_nChannels,lgd_cc_tune,&actualRun)){
    fprintf(stderr,"Error reading LGD tune constants for run %d!\n",event->runNo);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	  "lgdTuneLoad successful for run %d, real run=%d.\n",
	  event->runNo,actualRun);
#endif

  /* Load pedestals, calibration constants, thresholds for BGV  */

  if (bgv_ped)
    free(bgv_ped);
  bgv_ped = malloc(N_BGV_CHANNELS*sizeof(float));
  if (bgv_thresh)
    free(bgv_thresh);
  bgv_thresh = malloc(N_BGV_CHANNELS*sizeof(int));
  if (pedestalLoad(PED_BGV,event->runNo,N_BGV_CHANNELS,bgv_ped,bgv_thresh,&actualRun))
    return(1);

  if (bgv_cc)
    free(bgv_cc);
  bgv_cc = malloc(N_BGV_CHANNELS*sizeof(float));
  if (ccLoad(CC_BGV,event->runNo,N_BGV_CHANNELS,bgv_cc,&actualRun))
    return(1);

  if (bgv_tdc_offset)
    free(bgv_tdc_offset);
  bgv_tdc_offset = malloc(N_BGV_CHANNELS*sizeof(float));
  if (tdcOffsetLoad(TDC_OFFSET_BGV,event->runNo,N_BGV_CHANNELS,
		   bgv_tdc_offset,&bgv_t0_offset,&actualRun)){
    fprintf(stderr,
	    "Error reading BGV TDC offset run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "tdcOffsetLoad (BGV) success: run %d, real run=%d\n",
	    event->runNo,actualRun);
  }
#endif

  if (bgv_offset_ave)
    free(bgv_offset_ave);
  bgv_offset_ave = malloc(N_BGV_COUNTERS*sizeof(float));
  if (offsetAveLoadBGV(event->runNo,N_BGV_COUNTERS,bgv_offset_ave,
		      &bgv_t0_offset,&actualRun)){
    fprintf(stderr,
	    "Error reading BGV ave offset run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "aveOffsetLoadBGV success: run %d, real run=%d.\n",
	    event->runNo,actualRun);
  }
#endif

  if (bgv_mean_diff)
    free(bgv_mean_diff);
  bgv_mean_diff = malloc(N_BGV_COUNTERS*sizeof(float));
  if (meanDiffLoadBGV(event->runNo,N_BGV_COUNTERS,bgv_mean_diff,
		     &actualRun)){
    fprintf(stderr,
	    "Error reading BGV mean diff run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "meanDiffLoadBGV success: run %d, real run=%d.\n",
	    event->runNo,actualRun);
  }
#endif

  if (bgv_US_offset)
    free(bgv_US_offset);
  bgv_US_offset = malloc(N_BGV_COUNTERS*sizeof(float));
  if (offsetUSLoadBGV(event->runNo,N_BGV_COUNTERS,bgv_US_offset,
		     &bgv_t0_offset,&actualRun)){
    fprintf(stderr,
	    "Error reading BGV US offset run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "offsetUSLoadBGV success: run %d, real run=%d.\n",
	    event->runNo,actualRun);
  }
#endif

  if (bgv_DS_offset)
    free(bgv_DS_offset);
  bgv_DS_offset = malloc(N_BGV_COUNTERS*sizeof(float));
  if (offsetDSLoadBGV(event->runNo,N_BGV_COUNTERS,bgv_DS_offset,
		     &bgv_t0_offset,&actualRun)){
    fprintf(stderr,
	    "Error reading BGV DS offset run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "offsetDSLoadBGV success: run %d, real run=%d.\n",
	    event->runNo,actualRun);
  }
#endif

  /* Load calibration constants for tagger */

  if (tag_tdc_offset)
    free(tag_tdc_offset);
  tag_tdc_offset = malloc(N_TAG_CHANNELS*sizeof(float));
#if VERBOSE
  fprintf(stderr,"Tagger TDC offsets being read from map.\n");
#endif
  if (tdcOffsetLoad(TDC_OFFSET_TAG,event->runNo,N_TAG_CHANNELS,
		   tag_tdc_offset,&tag_t0_offset,&actualRun)){
    fprintf(stderr,
	    "Error reading TAG TDC offset run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "tdcOffsetLoad (TAG) success: run %d, real run=%d\n",
	    event->runNo,actualRun);
  }
#endif

  /* Load pedestals, calibration constants, thresholds for CPV */

  if (cpv_ped)
    free(cpv_ped);
  cpv_ped = malloc(N_CPV_CHANNELS*sizeof(float));
  if (cpv_thresh)
    free(cpv_thresh);
  cpv_thresh = malloc(N_CPV_CHANNELS*sizeof(int));
  if (pedestalLoad(PED_CPV,event->runNo,N_CPV_CHANNELS,cpv_ped,cpv_thresh,&actualRun))
    return(1);

  if (cpv_cc)
    free(cpv_cc);
  cpv_cc = malloc(N_CPV_CHANNELS*sizeof(float));
  if (ccLoad(CC_CPV,event->runNo,N_CPV_CHANNELS,cpv_cc,&actualRun))
    return(1);

  if (cpv_tdc_offset)
    free(cpv_tdc_offset);
  cpv_tdc_offset = malloc(N_CPV_CHANNELS*sizeof(float));
  if (tdcOffsetLoad(TDC_OFFSET_CPV,event->runNo,N_CPV_CHANNELS,
		   cpv_tdc_offset,&cpv_t0_offset,&actualRun)){
    fprintf(stderr,
	    "Error reading CPV TDC offset run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "tdcOffsetLoad (CPV) success: run %d, real run=%d\n",
	    event->runNo,actualRun);
  }
#endif

  /* Load pedestals, calibration constants, thresholds for UPV */

  if (upv_ped)
    free(upv_ped);
  upv_ped = malloc(N_UPV_CHANNELS*sizeof(float));
  if (upv_thresh)
    free(upv_thresh);
  upv_thresh = malloc(N_UPV_CHANNELS*sizeof(int));
  if (pedestalLoad(PED_UPV,event->runNo,N_UPV_CHANNELS,upv_ped,upv_thresh,&actualRun))
    return(1);

  if (upv_cc)
    free(upv_cc);
  upv_cc = malloc(N_UPV_CHANNELS*sizeof(float));
  if (ccLoad(CC_UPV,event->runNo,N_UPV_CHANNELS,upv_cc,&actualRun))
    return(1);

  if (upv_tdc_offset)
    free(upv_tdc_offset);
  upv_tdc_offset = malloc(N_UPV_CHANNELS*sizeof(float));
  if (tdcOffsetLoad(TDC_OFFSET_UPV,event->runNo,N_UPV_CHANNELS,
		   upv_tdc_offset,&upv_t0_offset,&actualRun)){
    fprintf(stderr,
	    "Error reading UPV TDC offset run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "tdcOffsetLoad (UPV) success: run %d, real run=%d\n",
	    event->runNo,actualRun);
  }
#endif
 
  /* Load constants needed for forming BSD pixels */

#if VERBOSE
  fprintf(stderr,"Loading pixel position information.\n");
#endif
  if (pixelPositionLoad(event->runNo,&bsd_ring_midz,&actualRun)){
    fprintf(stderr,
	      "Error reading pixel Z position information run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "pixelPositionLoad (bsd_ring_midz) load success: run %d, real run=%d\n",
	      event->runNo,actualRun);
  }
#endif

#if VERBOSE
  fprintf(stderr,"Loading pixel gain information.\n");
#endif
  if (pixelGainLoad(event->runNo,&bsd_pixel_gain,&actualRun)){
    fprintf(stderr,
	      "Error reading pixel gain information run %d, real run %d.\n",
	      event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	      "pixelGainLoad (bsd_pixel_gain) load success: run %d, real run=%d\n",
	      event->runNo,actualRun);
  }
#endif

  /* Load tagger energy settings information */

#if VERBOSE
  fprintf(stderr,"Tagger energy settings being read from map.\n");
#endif
  if (taggerSettingsLoad(event->runNo,&tagger_energy_settings,
                                      &tagger_fraction_by_channel,
                                      &tagger_CPleakage_by_channel,
                      &actualRun)){
    fprintf(stderr,
	    "Error reading TAG energy settings run %d, real run %d.\n",
	    event->runNo,actualRun);
    return(1);
  }
#if VERBOSE
  else {
    fprintf(stderr,
	    "taggedEnergyLoad (TAG) success: run %d, real run=%d\n",
	    event->runNo,actualRun);
  }
#endif

  /* Load channel status arrays */

  if (channel_status_tag){
    free(channel_status_tag);
    channel_status_tag=NULL;
  }
  channel_status_tag=malloc(N_TAG_CHANNELS*sizeof(unsigned char));
  if (channel_status_upv){
    free(channel_status_upv);
    channel_status_upv=NULL;
  }
  channel_status_upv=malloc(N_UPV_CHANNELS*sizeof(unsigned char));
  if (channel_status_bsd){
    free(channel_status_bsd);
    channel_status_bsd=NULL;
  }
  channel_status_bsd=malloc(N_BSD_CHANNELS*sizeof(unsigned char));
  if (channel_status_bgv){
    free(channel_status_bgv);
    channel_status_bgv=NULL;
  }
  channel_status_bgv=malloc(N_BGV_CHANNELS*sizeof(unsigned char));
  if (channel_status_cpv){
    free(channel_status_cpv);
    channel_status_cpv=NULL;
  }
  channel_status_cpv=malloc(N_CPV_CHANNELS*sizeof(unsigned char));
  if (channel_status_lgd){
    free(channel_status_lgd);
    channel_status_lgd=NULL;
  }
  channel_status_lgd=malloc(lgd_nChannels*sizeof(unsigned char));
  
    
  /* Setup the clusterizer */

  lgdSetup(3,0.15,0.05,LGDCLUS_LOG2,0);
  bgvSetup(0.02,0.005);           /* set the primary energy threshold to 20 MeV
                                     and the secondary threshold to 5 MeV */
  
  if (channelStatusLoad(event->runNo,
		       channel_status_tag,channel_status_upv,
		       channel_status_bsd,channel_status_bgv,
		       channel_status_cpv,channel_status_lgd,
		       actualRunList)){
    fprintf(stderr,"failed to load channel status maps!\n");
    return(1);
  }
  else {
#if VERBOSE
    fprintf(stderr,"Succesfully loaded all channel status maps\n");
    fprintf(stderr,"From runs: %4d %4d %4d %4d %4d %4d\n",
	    actualRunList[0],actualRunList[1],actualRunList[2],
	    actualRunList[3],actualRunList[4],actualRunList[5]);
#endif
#ifdef DEBUGGING
    fprintf(stderr,"LGD status map:\n");
    for(i=27;i>=0;i--){ /* i is row number */
      fprintf(stderr,"row %2d ",i);
      for(j=0;j<lgd_nCols;j++){ /* j is column number */
	ch=(i*lgd_nRows)+j;
	if (channel_status_lgd[ch] & CHANNEL_STATUS_IS_A_HOLE){
	  fprintf(stderr,"  ");
	}
	else if (channel_status_lgd[ch]){
	  fprintf(stderr," *");
	}
	else{
	  fprintf(stderr," .");
	}
      }
      fprintf(stderr,"\n");
    }
#endif
  }
  
  if (setup_lgd_cluster_cleanup(NULL)){
    fprintf(stderr,"failed to set up cluster cleanup.\n");
    return(1);
  }
  if (setup_makePhotons(event)){
    fprintf(stderr,"There was an error setting up make_photons.\n");
    return(1);
  }
  if (setup_makeMesons(event,NULL)){
    fprintf(stderr,"There was an error setting up make_mesons.\n");
    return(1);
  }

  /* Load a random number sequence for curing spikes in histograms */

  if (randomNo == NULL) {
    int i;
    
    randomNo = malloc(randomNos*sizeof(float));

    for(i=0; i<randomNos; i++)
      randomNo[i]=drand48();
  }

#if VERBOSE
  fprintf(stderr,"setup_makeHits successfully finished.\n");
#endif

  return(0);
}

#ifdef OLD_NONLINEAR_CORRECTION
/*---------------------------------------------------------------------------*/
/*  This is a routine to override the default values of the non-linearity    */
/*  correction parameters. Must be called some time after setup_makeHits     */
/*  in order for it to work properly. The function that calls this also      */
/*  needs to include makeHits.h     D.S. Armstrong                           */

void setupNonLinearity(float new_beta, float new_epsilon) {

  beta = new_beta;
  epsilon = new_epsilon;
  return;
}
#endif
