#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calibration.h>
#include <map_manager.h>
#include <detector_params.h>
#include <itypes.h>
#include <makeHits.h>
#include <lgdGeom.h>

/* #define DEBUGGING 1 */

static char *mapFile=NULL;
static char *map = "/maps/detector_params.map";

int meanDiffLoadBGV(int runNo, int nChannels, 
		    float *mean_diffs, int *actualRun)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"meanDiffLoadBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"mean_diff_%d",nChannels);
  if(map_get_float(mapFile,"BGV",str,
		   nChannels,mean_diffs,runNo,actualRun)){
    fprintf(stderr,"meanDiffLoadBGV: Failed to load BGV %s for run %d\n",
		    str,runNo);
    return(1);
  }
  return(0);
}

int offsetAveLoadBGV(int runNo, int nChannels, 
		     float *ave_offset, float *t0_offset, int *actualRun)
{
  char *dir;
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"offsetAveLoadBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"ave_offset_%d",nChannels);
  if(map_get_float(mapFile,"BGV",str,
		   nChannels,ave_offset,runNo,actualRun)){
    fprintf(stderr,"offsetAveLoadBGV: Failed to load BGV %s for run %d\n",
		    str,runNo);
    return(1);
  }
  sprintf(str,"t0_offset");
  if(map_get_float(mapFile,"BGV",str,
		   1,t0_offset,runNo,actualRun)){
    fprintf(stderr,"offsetAveLoadBGV: Failed to load BGV %s for run %d\n",
		    str,runNo);
    return(1);
  }
  return(0);
}

int offsetUSLoadBGV(int runNo, int nChannels,
		    float *US_offset, float *t0_offset, int *actualRun)
{
  char *dir;
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"offsetUSLoadBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"US_offset");
  if(map_get_float(mapFile,"BGV",str,
		   nChannels,US_offset,runNo,actualRun)){
    fprintf(stderr,"offset_USLoadBGV: Failed to load BGV %s for run %d\n",
		    str,runNo);
    return(1);
  }
  sprintf(str,"t0_offset");
  if(map_get_float(mapFile,"BGV",str,
		   1,t0_offset,runNo,actualRun)){
    fprintf(stderr,"offset_USLoadBGV: Failed to load BGV %s for run %d\n",
		    str,runNo);
    return(1);
  }
  return(0);
}

int offsetDSLoadBGV(int runNo, int nChannels,
		    float *DS_offset, float *t0_offset, int *actualRun)
{
  char *dir;
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"offsetDSLoadBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"DS_offset");
  if(map_get_float(mapFile,"BGV",str,
		   nChannels,DS_offset,runNo,actualRun)){
    fprintf(stderr,"offsetDSLoadBGV: Failed to load BGV %s for run %d\n",
	    str,runNo);
    return(1);
  }
  sprintf(str,"t0_offset");
  if(map_get_float(mapFile,"BGV",str,
		   1,t0_offset,runNo,actualRun)){
    fprintf(stderr,"offsetDSLoadBGV: Failed to load BGV %s for run %d\n",
	    str,runNo);
    return(1);
  }
  return(0);
}

int attenLengthLoadBGV(int runNo, int nChannels, 
		       float *atten_lengths, int *actualRun)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"attenLengthLoadBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"atten_length_%d",nChannels);
  if (map_get_float(mapFile,"BGV",str,
		   nChannels,atten_lengths,runNo,actualRun)){
    fprintf(stderr,"attenLengthLoadBGV: Failed to load BGV %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int attenLengthLoadBSD(int runNo, int nChannels, 
		       float *atten_lengths, int *actualRun)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"attenLengthLoadBSD: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"atten_length_%d",nChannels);
  if(map_get_float(mapFile,"BSD",str,
		   nChannels,atten_lengths,runNo,actualRun)){
    fprintf(stderr,"attenLengthLoadBSD: Failed to load BSD %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int tdcBoundLoad(int detector, int runNo, int nChannels, 
		 int *bound, int *actualRun)
{
  char *dir;
  char str[20];
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"tdcBoundLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  switch(detector){
  case TDC_UPPER_BOUND_BSD:
    sprintf(str,"high_edge_%d",nChannels);
    if(map_get_int(mapFile,"BSD",str,nChannels,bound,runNo,actualRun)){
      fprintf(stderr,"tdcBoundLoad: load for BSD %s, run %d, failed\n",str,runNo);
      return(1);
    }
    break;
  case TDC_LOWER_BOUND_BSD:
    sprintf(str,"low_edge_%d",nChannels);
    if(map_get_int(mapFile,"BSD",str,nChannels,bound,runNo,actualRun)){
      fprintf(stderr,"tdcBoundLoad: load for BSD %s, run %d, failed\n",str,runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"tdcBoundLoad: Unknown detector (0x%x)\n",detector);
    return(1);
  }
  return(0);
}

int pixelTimingLoad(int runNo, 
		    float *pixel_LR_offset, float *pixel_LR_width,
		    float *pixel_LS_offset, float *pixel_LS_width,
		    float *pixel_RS_offset, float *pixel_RS_width,
		    int *actualRun)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pixelTimingLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  
  sprintf(str,"pixel_LR_diff");
  if(map_get_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LR_offset,runNo,actualRun)){
    fprintf(stderr,"pixelTimingLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_LR_width");
  if(map_get_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LR_width,runNo,actualRun)){
    fprintf(stderr,"pixelTimingLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_LS_diff");
  if(map_get_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LS_offset,runNo,actualRun)){
    fprintf(stderr,"pixelTimingLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_LS_width");
  if(map_get_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LS_width,runNo,actualRun)){
    fprintf(stderr,"pixelTimingLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_RS_diff");
  if(map_get_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_RS_offset,runNo,actualRun)){
    fprintf(stderr,"pixelTimingLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_RS_width");
  if(map_get_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_RS_width,runNo,actualRun)){
    fprintf(stderr,"pixelTimingLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int pixelPositionLoad(int runNo, float *pixel_Z_by_ring, int *actualRun)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pixelPositionLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  
  sprintf(str,"pixel_Z_by_ring");
  if(map_get_float(mapFile,"BSD",str,BSD_N_RINGS,
		   pixel_Z_by_ring,runNo,actualRun)){
    fprintf(stderr,"pixelTimingLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int pixelGainLoad(int runNo, float *pixel_gain, int *actualRun)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pixelGainLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  
  sprintf(str,"pixel_gain_384");
  if(map_get_float(mapFile,"BSD",str,384,
		   pixel_gain,runNo,actualRun)){
    fprintf(stderr,"pixelGainLoad: Failed to load BSD %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int taggerSettingsLoad(int runNo, float *tagger_settings,
                       float *fraction_by_channel,
                       float *CPleakage_by_channel,
                       int *actualRun)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"taggerSettingsLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  
  sprintf(str,"taggedEnergy");
  if(map_get_float(mapFile,"TAG",str,3,
		   tagger_settings,runNo,actualRun)){
    fprintf(stderr,"taggerSettingsLoad: Failed to load TAG %s for run %d.\n",str,runNo);
    return(1);
  }

  sprintf(str,"taggedFraction");
  if(map_get_float(mapFile,"TAG",str,N_TAGGER_ENERGY_CHANNELS,
		   fraction_by_channel,runNo,actualRun)){
    fprintf(stderr,"taggerSettingsLoad: Failed to load TAG %s for run %d.\n",str,runNo);
    return(1);
  }

  sprintf(str,"taggedCPleakage");
  if(map_get_float(mapFile,"TAG",str,N_TAGGER_ENERGY_CHANNELS,
		   CPleakage_by_channel,runNo,actualRun)){
    fprintf(stderr,"taggerSettingsLoad: Failed to load TAG %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int tdcOffsetLoad(int detector, int runNo, int nChannels,
	       	  float *tdc_offset, float *t0_offset, int *actualRun)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"tdcOffsetLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  switch(detector){
  case TDC_OFFSET_TAG:
    sprintf(str,"tdc_offset_%d",nChannels);
    if(map_get_float(mapFile,"TAG",str,nChannels,tdc_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load TAG %s for run %d\n",str,runNo);
      return(1);
    }
    sprintf(str,"t0_offset");
    if(map_get_float(mapFile,"TAG",str,1,t0_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load TAG %s for run %d\n",str,runNo);
      return(1);
    }
    break;
  case TDC_OFFSET_BGV:
    sprintf(str,"tdc_offset_%d",nChannels);
    if(map_get_float(mapFile,"BGV",str,nChannels,tdc_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load BGV %s for run %d\n",str,runNo);
      return(1);
    }
    sprintf(str,"t0_offset");
    if(map_get_float(mapFile,"BGV",str,1,t0_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load BGV %s for run %d\n",str,runNo);
      return(1);
    }
    break;
  case TDC_OFFSET_BSD:
    sprintf(str,"tdc_offset_%d",nChannels);
    if(map_get_float(mapFile,"BSD",str,nChannels,tdc_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load BSD %s for run %d\n",str,runNo);
      return(1);
    }
    sprintf(str,"t0_offset");
    if(map_get_float(mapFile,"BSD",str,1,t0_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load BSD %s for run %d\n",str,runNo);
      return(1);
    }
    break;
  case TDC_OFFSET_CPV:
    sprintf(str,"tdc_offset_%d",nChannels);
    if(map_get_float(mapFile,"CPV",str,nChannels,tdc_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load CPV %s for run %d\n",str,runNo);
      return(1);
    }
    sprintf(str,"t0_offset");
    if(map_get_float(mapFile,"CPV",str,1,t0_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load CPV %s for run %d\n",str,runNo);
      return(1);
    }
    break;
  case TDC_OFFSET_UPV:
    sprintf(str,"tdc_offset_%d",nChannels);
    if(map_get_float(mapFile,"UPV",str,nChannels,tdc_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load UPV %s for run %d\n",str,runNo);
      return(1);
    }
    sprintf(str,"t0_offset");
    if(map_get_float(mapFile,"UPV",str,1,t0_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load UPV %s for run %d\n",str,runNo);
      return(1);
    }
    break;
  case TIME_OFFSET_BSD_TAG:
    sprintf(str,"bsd_tag_diff");
    if(map_get_float(mapFile,"TAG",str,nChannels,tdc_offset,runNo,actualRun)){
      fprintf(stderr,"tdcOffsetLoad: Failed to load %s for run %d\n",str,runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"tdcOffsetLoad : Unknown detector (0x%x)\n",detector);
    return(1);
  }  
  return(0);
}

int channelStatusLoad(int runNo, 
		      unsigned char *channel_status_tag,
		      unsigned char *channel_status_upv,
		      unsigned char *channel_status_bsd,
		      unsigned char *channel_status_bgv,
		      unsigned char *channel_status_cpv,
		      unsigned char *channel_status_lgd,
		      int *actual_run_list)
{
  char *dir; 
  char str[20];
  int olength;
  int actualRun;
  lgdGeom_t geometry_query;
  int i;
  int return_val;

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"channelStatusLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"channel_status");

  if(map_get_char(mapFile,"TAG",str,N_TAG_CHANNELS,
		  (char*)channel_status_tag,runNo,actual_run_list+0,&olength)){
    fprintf(stderr,"channelStatusLoad: Failed to load channel_status_tag for run %d\n",runNo);
    return(1);
  }
  if(map_get_char(mapFile,"UPV",str,N_UPV_CHANNELS,
		  (char*)channel_status_upv,runNo,actual_run_list+1,&olength)){
    fprintf(stderr,"channelStatusLoad: Failed to load channel_status_upv for run %d\n",runNo);
    return(1);
  }
  if(map_get_char(mapFile,"BSD",str,N_BSD_CHANNELS,
		  (char*)channel_status_bsd,runNo,actual_run_list+2,&olength)){
    fprintf(stderr,"channelStatusLoad: Failed to load channel_status_bsd for run %d\n",runNo);
    return(1);
  }
  if(map_get_char(mapFile,"BGV",str,N_BGV_CHANNELS,
		  (char*)channel_status_bgv,runNo,actual_run_list+3,&olength)){
    fprintf(stderr,"channelStatusLoad: Failed to load channel_status_bgv for run %d\n",runNo);
    return(1);
  }
  if(map_get_char(mapFile,"CPV",str,N_CPV_CHANNELS,
		  (char*)channel_status_cpv,runNo,actual_run_list+4,&olength)){
    fprintf(stderr,"channelStatusLoad: Failed to load channel_status_cpv for run %d\n",runNo);
    return(1);
  }
  if(map_get_char(mapFile,"LGD",str,lgd_nChannels,
		  (char*)channel_status_lgd,runNo,actual_run_list+5,&olength)){
    fprintf(stderr,"channelStatusLoad: Failed to load channel_status_lgd for run %d\n",runNo);
    return(1);
  }
  /* now populates LGD holes into array using lgd Geometry libraries */
  for(i=0;i<lgd_nChannels;i++){
    geometry_query.channel=i;
    return_val=lgdGetGeom(0,&geometry_query,LGDGEOM_BYCHANNEL);
    switch(return_val){
    case LGDGEOM_OK:
      continue;
      break;
    case LGDGEOM_ISAHOLE:
      (channel_status_lgd[i])|=CHANNEL_STATUS_IS_A_HOLE;
#ifdef DEBUGGING 
      fprintf(stderr,"lgd channel %3d is a hole\n",i);
#endif
      continue;
      break;
    default:
      fprintf(stderr,"channelStatusLoad: lgdGetGeom returned unexpected result:\n");
      fprintf(stderr,"channel %d, result=%d\n",
	      i,return_val);
      return(1);
    }
  } /* for i over all LGD channels for filling hole array */
  return(0);
}

int meanDiffSaveBGV(int runNo, int nChannels, float *mean_diffs)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"meanDiffSaveBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  sprintf(str,"mean_diff_%d",nChannels);
  if(map_put_float(mapFile,"BGV",str,
		   nChannels,mean_diffs,runNo)){
    fprintf(stderr,
	    "meanDiffSaveBGV: Failed to upload BGV %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int offsetAveSaveBGV(int runNo, int nChannels,
	             float *ave_offset, float *t0_offset)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"offsetAveSaveBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  sprintf(str,"ave_offset_%d",nChannels);
  if(map_put_float(mapFile,"BGV",str,
		   nChannels,ave_offset,runNo)){
    fprintf(stderr,
	    "offsetAveSaveBGV: Failed to upload BGV %s for run %d\n",str,runNo);
    return(1);
  }
  sprintf(str,"t0_offset");
  if(map_put_float(mapFile,"BGV",str,
		   1,t0_offset,runNo)){
    fprintf(stderr,"offsetAveSaveBGV: Failed to upload BGV %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int offsetUSSaveBGV(int runNo, int nChannels,
                    float *US_offset, float *t0_offset)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"offsetUSSaveBSD: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  sprintf(str,"US_offset");
  if(map_put_float(mapFile,"BGV",str,
		   nChannels,US_offset,runNo)){
    fprintf(stderr,"offsetUSSaveBGV: Failed to upload BGV %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int offsetDSSaveBGV(int runNo, int nChannels,
		    float *DS_offset, float *t0_offset)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"offsetDSSaveBSD: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  sprintf(str,"DS_offset");
  if(map_put_float(mapFile,"BGV",str,
		   nChannels,DS_offset,runNo)){
    fprintf(stderr,"offsetDSSaveBGV: Failed to upload BGV %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int attenLengthSaveBGV(int runNo, int nChannels, float *atten_length)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"attenLengthSaveBGV: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  sprintf(str,"atten_length_%d",nChannels);
  if(map_put_float(mapFile,"BGV",str,
		   nChannels,atten_length,runNo)){
    fprintf(stderr,"attenLengthSaveBGV: Failed to upload BGV %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int attenLengthSaveBSD(int runNo, int nChannels, float *atten_length)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"attenLengthSaveBSD: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  sprintf(str,"atten_length_%d",nChannels);
  if(map_put_float(mapFile,"BSD",str,
		   nChannels,atten_length,runNo)){
    fprintf(stderr,"attenLengthSaveBSD: Failed to upload BSD %s for run %d\n",str,runNo);
    return(1);
  }
  return(0);
}

int pixelTimingSave(int runNo, 
		    float *pixel_LR_offset, float *pixel_LR_width,
		    float *pixel_LS_offset, float *pixel_LS_width,
		    float *pixel_RS_offset, float *pixel_RS_width)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pixelTimingSave: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  
  sprintf(str,"pixel_LR_diff");
  if(map_put_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LR_offset,runNo)){
    fprintf(stderr,"pixelTimingSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_LR_width");
  if(map_put_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LR_width,runNo)){
    fprintf(stderr,"pixelTimingSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_LS_diff");
  if(map_put_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LS_offset,runNo)){
    fprintf(stderr,"pixelTimingSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_LS_width");
  if(map_put_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_LS_width,runNo)){
    fprintf(stderr,"pixelTimingSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_RS_diff");
  if(map_put_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_RS_offset,runNo)){
    fprintf(stderr,"pixelTimingSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
		    
  sprintf(str,"pixel_RS_width");
  if(map_put_float(mapFile,"BSD",str,N_DEFINED_BSD_PIXELS,
		   pixel_RS_width,runNo)){
    fprintf(stderr,"pixelTimingSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int pixelPositionSave(int runNo, float *pixel_Z_by_ring)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pixelPositionLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  
  sprintf(str,"pixel_Z_by_ring");
  if(map_put_float(mapFile,"BSD",str,BSD_N_RINGS,
		   pixel_Z_by_ring,runNo)){
    fprintf(stderr,"pixelTimingSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int pixelGainSave(int runNo, float *pixel_gain)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pixelGainLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  
  sprintf(str,"pixel_gain_384");
  if(map_put_float(mapFile,"BSD",str,384,
		   pixel_gain,runNo)){
    fprintf(stderr,"pixelGainSave: Failed to upload BSD %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int taggerSettingsSave(int runNo, float *tagger_settings,
                       float *fraction_by_channel)
{
  char *dir; 
  char str[20];

  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"taggerEnergyLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  
  sprintf(str,"taggedEnergy");
  if(map_put_float(mapFile,"TAG",str,3,
		   tagger_settings,runNo)){
    fprintf(stderr,"taggerEnergySave: Failed to upload TAG %s for run %d.\n",str,runNo);
    return(1);
  }

  sprintf(str,"taggedFraction");
  if(map_put_float(mapFile,"TAG",str,N_TAGGER_ENERGY_CHANNELS,
		   fraction_by_channel,runNo)){
    fprintf(stderr,"taggerEnergySave: Failed to upload TAG %s for run %d.\n",str,runNo);
    return(1);
  }
  return(0);
}

int channelStatusSave(int detector, int runNo,
		  int nChannels, unsigned char *channel_status)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"channelStatusSave: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  sprintf(str,"channel_status");        

  switch(detector){
    
  case CHANNEL_STATUS_TAG:
    if(N_TAG_CHANNELS!=nChannels){
      fprintf(stderr,"channelStatusSave(TAG): nChannels wrong!\n");
      fprintf(stderr,"nChannels=%d, should be %d\n",
	      nChannels,N_TAG_CHANNELS);
      return(1);
    }
    if(map_put_char(mapFile,"TAG",str,nChannels,(char*)channel_status
		    ,runNo)){
      fprintf(stderr,"channelStatusSave: Failed to save channel_status_tag for run %d\n",runNo);
      return(1);
    }
    break;
    
  case CHANNEL_STATUS_UPV:
    if(N_UPV_CHANNELS!=nChannels){
      fprintf(stderr,"channelStatusSave(UPV): nChannels wrong!\n");
      fprintf(stderr,"nChannels=%d, should be %d\n",
	      nChannels,N_UPV_CHANNELS);
      return(1);
    }
    if(map_put_char(mapFile,"UPV",str,nChannels,(char*)channel_status
		    ,runNo)){
      fprintf(stderr,"channelStatusSave: Failed to save channel_status_upv for run %d\n",runNo);
      return(1);
    }
    break;
    
  case CHANNEL_STATUS_BSD:
    if(N_BSD_CHANNELS!=nChannels){
      fprintf(stderr,"channelStatusSave(BSD): nChannels wrong!\n");
      fprintf(stderr,"nChannels=%d, should be %d\n",
	      nChannels,N_BSD_CHANNELS);
      return(1);
    }
    if(map_put_char(mapFile,"BSD",str,nChannels,(char*)channel_status
		    ,runNo)){
      fprintf(stderr,"channelStatusSave: Failed to save channel_status_bsd for run %d\n",runNo);
      return(1);
    }
    break;
    
  case CHANNEL_STATUS_BGV:
    if(N_BGV_CHANNELS!=nChannels){
      fprintf(stderr,"channelStatusSave(BGV): nChannels wrong!\n");
      fprintf(stderr,"nChannels=%d, should be %d\n",
	      nChannels,N_BGV_CHANNELS);
      return(1);
    }
    if(map_put_char(mapFile,"BGV",str,nChannels,(char*)channel_status
		    ,runNo)){
      fprintf(stderr,"channelStatusSave: Failed to save channel_status_bgv for run %d\n",runNo);
      return(1);
    }
    break;
    
  case CHANNEL_STATUS_CPV:
    if(N_CPV_CHANNELS!=nChannels){
      fprintf(stderr,"channelStatusSave(CPV): nChannels wrong!\n");
      fprintf(stderr,"nChannels=%d, should be %d\n",
	      nChannels,N_CPV_CHANNELS);
      return(1);
    }
    if(map_put_char(mapFile,"CPV",str,nChannels,(char*)channel_status
		    ,runNo)){
      fprintf(stderr,"channelStatusSave: Failed to save channel_status_cpv for run %d\n",runNo);
      return(1);
    }
    break;
    
  case CHANNEL_STATUS_LGD:
    if(lgd_nChannels!=nChannels){
      fprintf(stderr,"channelStatusSave(LGD): nChannels wrong!\n");
      fprintf(stderr,"nChannels=%d, should be %d\n",
	      nChannels,lgd_nChannels);
      return(1);
    }
    if(map_put_char(mapFile,"LGD",str,nChannels,(char*)channel_status
		    ,runNo)){
      fprintf(stderr,"channelStatusSave: Failed to save channel_status_lgd for run %d\n",runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"channelStatusSave: ERROR!\n");
    fprintf(stderr,"Tried to save channel status for detector %d\n",
	    detector);
    fprintf(stderr,"There is no such thing!\n");
    return(1);
  }
  fprintf(stderr,"Successfully saved channel_status array.\n");
  return(0);
}
	    

int tdcOffsetSave(int detector, int runNo, 
		  int nChannels, float *tdc_offset, float *t0_offset)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,
	      "tdcOffsetSave: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  switch(detector){
  case TDC_OFFSET_BSD:
    sprintf(str,"tdc_offset_%d",nChannels);
    if(map_put_float(mapFile,"BSD",str,
		     nChannels,tdc_offset,runNo)){
      fprintf(stderr,"tdcOffsetSave: Failed save BSD %s run %d.\n",str,runNo);
      return(1);
    }
    sprintf(str,"t0_offset");
    if(map_put_float(mapFile,"BSD",str,
		     1,t0_offset,runNo)){
      fprintf(stderr,"tdcOffsetSave: Failed save BSD %s run %d.\n",str,runNo);
      return(1);
    }
    break;
  case TDC_OFFSET_TAG:
    sprintf(str,"tdc_offset_%d",nChannels);
    if(map_put_float(mapFile,"TAG",str,
		     nChannels,tdc_offset,runNo)){
      fprintf(stderr,"tdcOffsetSave: Failed save TAG %s run %d.\n",str,runNo);
      return(1);
    }
    sprintf(str,"t0_offset");
    if(map_put_float(mapFile,"TAG",str,
		     1,t0_offset,runNo)){
      fprintf(stderr,"tdcOffsetSave: Failed save TAG %s run %d.\n",str,runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"tdcOffsetSave : Unknown detector (0x%x)\n",detector);
    return(1);
  }  
  return(0);
}

int tdcBoundSave(int detector, int runNo, int nChannels, int *bound)
{
  char *dir;
  char str[20];
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"tdcBoundSave: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  switch(detector){
  case TDC_UPPER_BOUND_BSD:
    sprintf(str,"high_edge_%d",nChannels);
    if(map_put_int(mapFile,"BSD",str,nChannels,bound,runNo)){
      fprintf(stderr,"tdcBoundSave: upload of BSD %s, run %d, failed\n",str,runNo);
      return(1);
    }
    break;
  case TDC_LOWER_BOUND_BSD:
    sprintf(str,"low_edge_%d",nChannels);
    if(map_put_int(mapFile,"BSD",str,nChannels,bound,runNo)){
      fprintf(stderr,"tdcBoundSave: upload of BSD %s, run %d, failed\n",str,runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"tdcBoundSave: Unknown detector (0x%x)\n",detector);
    return(1);
  }
  return(0);
}
