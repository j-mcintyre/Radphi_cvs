/*
  $Log: pedestals.c,v $
  Revision 1.8  2004/09/30 19:26:16  radphi
  detector_params.c, pedestals.c
     - modified to pull constants from the maps that were formerly
       hard-wired in header files -rtj-

  Revision 1.7  2001/04/04 20:22:32  radphi
  Changes made by crsteffe@mantrid00
  has provisoins for channel map vetoing

  Revision 1.6  2001/03/06 15:59:15  radphi
  Changes made by crsteffe@mantrid00
  updated to use new threshold algorithm

  Revision 1.5  1999/07/04 18:31:31  radphi
  Changes made by gevmage@jlabs1
  ooops...

  Revision 1.4  1999/07/04 18:17:58  radphi
  Changes made by gevmage@jlabs1
  I've added the 1999 detectors, as well as the ADC control code

  Revision 1.3  1997/06/02 20:40:54  radphi
  Changes made by lfcrob@dustbunny
  Added rcsid[]

 * Revision 1.2  1997/05/29  20:37:02  radphi
 * Changes made by lfcrob@dustbunny
 * Added support for ESUM,RPD, and CPV
 *
 * Revision 1.1  1997/05/16  19:19:16  radphi
 * Initial revision by lfcrob@jlabs2
 * First stab pedestal calculator
 * Only does LGD right now
 *
  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pedestals.h>
#include <map_manager.h>
#include <itypes.h>

static const char rcsid[]="$Id: pedestals.c,v 1.8 2004/09/30 19:26:16 radphi Exp $";

#define SIMPLIFIED_THRESHOLDS 1

static char *mapFile=NULL;
static char *mapFile_adc_control=NULL;
static char *map = "/maps/pedestals.map";
static char *map_adc_control = "/maps/adc_control.map";

int pedestalLoad(int detector, int runNo, int nChannels, float *pedestal, int *threshold, int *actualRun)
{
  char *dir; 
  char str[20];
  int i;
  int margin;
  int actualRun2;
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pedestalLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  switch(detector){
  case PED_LGD:
    sprintf(str,"ped_%d",nChannels);        
    if(map_get_float(mapFile,"LGD",str,nChannels,pedestal,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load LGD pedestals for run %d\n",runNo);
      return(1);
    }
#ifndef SIMPLIFIED_THRESHOLDS
    sprintf(str,"thresh_%d",nChannels);        
    if(map_get_int(mapFile,"LGD",str,nChannels,threshold,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load LGD thresholds for run %d\n",runNo);
      return(1);
    }
#endif
#ifdef SIMPLIFIED_THRESHOLDS
    sprintf(str,"LGD");
#endif
    break;
  case PED_UPV:
    sprintf(str,"ped_%d",nChannels);        
    if(map_get_float(mapFile,"UPV",str,nChannels,pedestal,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load UPV pedestals for run %d\n",runNo);
      return(1);
    }
#ifndef SIMPLIFIED_THRESHOLDS
    sprintf(str,"thresh_%d",nChannels);        
    if(map_get_int(mapFile,"UPV",str,nChannels,threshold,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load UPV thresholds for run %d\n",runNo);
      return(1);
    }
#endif
#ifdef SIMPLIFIED_THRESHOLDS
    sprintf(str,"UPV");
#endif
    break;
  case PED_BSD:
    sprintf(str,"ped_%d",nChannels);        
    if(map_get_float(mapFile,"BSD",str,nChannels,pedestal,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load BSD pedestals for run %d\n",runNo);
      return(1);
    }
#ifndef SIMPLIFIED_THRESHOLDS
    sprintf(str,"thresh_%d",nChannels);        
    if(map_get_int(mapFile,"BSD",str,nChannels,threshold,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load BSD thresholds for run %d\n",runNo);
      return(1);
    }
#endif
#ifdef SIMPLIFIED_THRESHOLDS
    sprintf(str,"BSD");
#endif
    break;
  case PED_BGV:
    sprintf(str,"ped_%d",nChannels);        
    if(map_get_float(mapFile,"BGV",str,nChannels,pedestal,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load BGV pedestals for run %d\n",runNo);
      return(1);
    }
#ifndef SIMPLIFIED_THRESHOLDS
   sprintf(str,"thresh_%d",nChannels);        
    if(map_get_int(mapFile,"BGV",str,nChannels,threshold,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load BGV thresholds for run %d\n",runNo);
      return(1);
    }
#endif
#ifdef SIMPLIFIED_THRESHOLDS
    sprintf(str,"BGV");
#endif
    break;
  case PED_CPV:
    sprintf(str,"ped_%d",nChannels);        
    if(map_get_float(mapFile,"CPV",str,nChannels,pedestal,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load CPV pedestals for run %d\n",runNo);
      return(1);
    }
#ifndef SIMPLIFIED_THRESHOLDS
    sprintf(str,"thresh_%d",nChannels);        
    if(map_get_int(mapFile,"CPV",str,nChannels,threshold,runNo,actualRun)){
      fprintf(stderr,"pedestalLoad: Failed to load CPV thresholds for run %d\n",runNo);
      return(1);
    }
#endif
#ifdef SIMPLIFIED_THRESHOLDS
    sprintf(str,"CPV");
#endif
    break;
  default:
    fprintf(stderr,"pedestalLoad : Unknown detector (0x%x)\n",detector);
    return(1);
  }
  if(map_get_int(mapFile,str,"thresh_margin",1,&margin,
		 runNo,&actualRun2)){
    fprintf(stderr,"pedestalLoad: Failed to get threshold margin!!!\n");
    return(1);
  }
#if VERBOSE
  fprintf(stderr,
	  "Utilizing new threshold algorithm; thresh = ped + constant\n");
#endif
  for(i=0;i<nChannels;i++){
    threshold[i]=(((int)(pedestal[i]))+
		  margin);
  }
  return(0);
}


int pedestalSave(int detector, int runNo, int nChannels, float *pedestal, int *threshold)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pedestalSave: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  if(!mapFile_adc_control){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"pedestalSave: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile_adc_control = malloc(strlen(dir)+strlen(map_adc_control)+1);
    sprintf(mapFile_adc_control,"%s%s",dir,map_adc_control);
  }
  runNo = abs(runNo);
  switch(detector){
  case PED_LGD:
    sprintf(str,"ped_%d",nChannels);        
    if(map_put_float(mapFile,"LGD",str,nChannels,pedestal,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save LGD pedestals for run %d\n",runNo);
      return(1);
    }
    sprintf(str,"thresh_%d",nChannels);
    if(map_put_int(mapFile,"LGD",str,nChannels,threshold,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save LGD thresholds for run %d\n",runNo);
      return(1);
    }
    break;

  case PED_UPV:
    sprintf(str,"ped_%d",nChannels);
    if(map_put_float(mapFile,"UPV",str,nChannels,pedestal,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save UPV pedestals for run %d\n",runNo);
      return(1);
    }
    sprintf(str,"thresh_%d",nChannels);
    if(map_put_int(mapFile,"UPV",str,nChannels,threshold,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save UPV thresholds for run %d\n",runNo);
      return(1);
    }
    break;

  case PED_BSD:
    sprintf(str,"ped_%d",nChannels);
    if(map_put_float(mapFile,"BSD",str,nChannels,pedestal,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save BSD pedestals for run %d\n",runNo);
      return(1);
    }
    sprintf(str,"thresh_%d",nChannels);
    if(map_put_int(mapFile,"BSD",str,nChannels,threshold,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save BSD thresholds for run %d\n",runNo);
      return(1);
    }
    break;

  case PED_BGV:
    sprintf(str,"ped_%d",nChannels);        
    if(map_put_float(mapFile,"BGV",str,nChannels,pedestal,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save BGV pedestals for run %d\n",runNo);
      return(1);
    }
    sprintf(str,"thresh_%d",nChannels);
    if(map_put_int(mapFile,"BGV",str,nChannels,threshold,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save BGV thresholds for run %d\n",runNo);
      return(1);
    }
    break;

  case PED_CPV:
    sprintf(str,"ped_%d",nChannels);        
    if(map_put_float(mapFile,"CPV",str,nChannels,pedestal,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save CPV pedestals for run %d\n",runNo);
      return(1);
    }
    sprintf(str,"thresh_%d",nChannels);
    if(map_put_int(mapFile,"CPV",str,nChannels,threshold,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save CPV thresholds for run %d\n",runNo);
      return(1);
    }
    break;

  case PED_ADCS_TOP:
    sprintf(str,"top_thresh",nChannels);
    if(map_put_int(mapFile_adc_control,"ADC_THRESH",str,nChannels,threshold,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save top crate ADC thresholds for run %d\n",runNo);
      return(1);
    }
    break;

  case PED_ADCS_BOTTOM:
    sprintf(str,"bottom_thresh",nChannels);
    if(map_put_int(mapFile_adc_control,"ADC_THRESH",str,nChannels,threshold,runNo)){
      fprintf(stderr,"pedestalSave: Failed to save bottom crate ADC thresholds for run %d\n",runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"pedestalSave : Unknown detector (0x%x)\n",detector);
    return(1);
  }  
  return(0);
}


