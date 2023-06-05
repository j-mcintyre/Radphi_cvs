/*
  $Log: detector_params.h,v $
  Revision 1.10  2003/01/14 21:08:53  radphi
  -detector_params.h : added t0 offsets for each detector class
  -makeHits.h : added t0 offset arrays for each detector type
  -rtj-

  Revision 1.9  2002/09/17 17:16:04  radphi
  - extended Craig's list of codes for LGD channel status to distinguish
    GOOD (0) from HOLE (1) and COLD (2) and HOT (4) blocks.  -rtj-

  Revision 1.8  2001/07/16 21:05:35  radphi
  Changes made by crsteffe@mantrid00
  updating ot make sure its curent

  Revision 1.7  2001/07/02 21:00:01  radphi
  Changes made by crsteffe@mantrid00
  updates that include stuff for new tagger clusters and recoil structures

  Revision 1.6  2001/03/16 20:59:05  radphi
  Changes made by crsteffe@mantrid00
  includes stuff for channel_status_* arrays

  Revision 1.5  2001/03/08 20:22:44  radphi
  Changes made by crsteffe@mantrid00
  contains definitoins for CHANNELSTATUS bit masks, first revision

  Revision 1.4  2001/01/24 16:44:18  radphi
  Changes made by crsteffe@mantrid00
  contains constants adn stuctures for recoil reconstruction

  Revision 1.3  2000/10/20 12:50:41  radphi
  Changes made by crsteffe@monolith.physics.indiana.edu
  current version from IU

  Revision 1.2  2000/06/28 19:17:45  radphi
  Changes made by radphi@ifarms2
  Now has all the function protyotypes

  Revision 1.1  2000/06/27 15:24:36  radphi
  Initial revision by radphi@urs3
  governs the detector_params.map software

  Revision 1.2  2000/06/18 06:37:07  radphi
  Changes made by jonesrt@grendl
  updates to ccLoad for bsd,bgv,upv groups -rtj

 * Revision 1.1  1997/06/02  18:36:59  radphi
 * Initial revision by lfcrob@dustbunny
 * Calibration constant routine headers
 *
  */

#ifndef DETECTOR_PARAMS_H_INCLUDED
#define DETECTOR_PARAMS_H_INCLUDED

#define TDC_OFFSET_RPD 0x2
#define TDC_OFFSET_CPV 0x4
#define TDC_OFFSET_UPV 0x10
#define TDC_OFFSET_BGV 0x20
#define TDC_OFFSET_BSD 0x40
#define TIME_OFFSET_BSD 0x41
#define TDC_OFFSET_TAG 0x80
#define TIME_OFFSET_BSD_TAG 0x81

#define TDC_UPPER_BOUND_BSD 0x40
#define TDC_LOWER_BOUND_BSD 0x41

#define N_DEFINED_BSD_PIXELS 192

#define CHANNEL_STATUS_TAG 0x01
#define CHANNEL_STATUS_UPV 0x02
#define CHANNEL_STATUS_BSD 0x04
#define CHANNEL_STATUS_BGV 0x08
#define CHANNEL_STATUS_CPV 0x10
#define CHANNEL_STATUS_LGD 0x20

/* these constants are used in channel status maps */
/* These values are going to be mapped into unsigned chars, */
/* so besides "GOOD", there should be up to 8 entries, each */
/* a hex number with only 1 bit equalling 1 */
#define CHANNEL_STATUS_GOOD          0x00
#define CHANNEL_STATUS_IS_A_HOLE     0x01
#define CHANNEL_STATUS_COLD          0x02
#define CHANNEL_STATUS_HOT           0x04

int meanDiffLoadBGV(int runNo, int nChannels, 
		    float *mean_diffs, int *actualRun);
int offsetAveLoadBGV(int runNo, int nChannels,
		     float *ave_offset, float *t0_offset, int *actualRun);
int attenLengthLoadBSD(int runNo, int nChannels, 
		       float *atten_lengths, int *actualRun);
int tdcOffsetLoad(int detector, int runNo, int nChannels, 
		  float *tdc_offset, float *t0_offset, int *actualRun);
int meanDiffSaveBGV(int runNo, int nChannels, float *mean_diffs);
int offsetAveSaveBGV(int runNo, int nChannels, float *ave_offset, float *t0_offset);
int attenLengthSaveBSD(int runNo, int nChannels, float *atten_length);
int tdcOffsetSave(int detector, int runNo, 
		  int nChannels, float *tdc_offset, float *t0_offset);
int tdcBoundLoad(int detector, int runNo, int nChannels, int *bound, int *actualRun);
int tdcBoundSave(int detector, int runNo, int nChannels, int *bound);
int pixelTimingLoad(int runNo, 
		    float *pixel_LR_offset, float *pixel_LR_width,
		    float *pixel_LS_offset, float *pixel_LS_width,
		    float *pixel_RS_offset, float *pixel_RS_width,
		    int *actualRun);
int channelStatusLoad(int runNo, 
		      unsigned char *channel_status_tag,
		      unsigned char *channel_status_upv,
		      unsigned char *channel_status_bsd,
		      unsigned char *channel_status_bgv,
		      unsigned char *channel_status_cpv,
		      unsigned char *channel_status_lgd,
		      int *real_run_list);
int channelStatusSave(int detector, int runNo,
		      int nChannels, unsigned char *channel_status);


#endif 
