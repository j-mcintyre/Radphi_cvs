/*

  $Log: calibration.c,v $
  Revision 1.7  2001/11/10 16:09:03  radphi
  - updated detector_params.c to load tdc offsets tables for cpv and upv
  - updated calibration.c to load cc tables for bgv, cpv and upv
  -rtj-

  Revision 1.6  2001/06/03 17:55:27  radphi
  Changes made by crsteffe@mantrid00
  has provisions for lgdtune constants

  Revision 1.5  2000/07/31 07:43:56  radphi
  Changes made by gevmage@urs3
  removed the abs from tdcOffsetSave

  Revision 1.4  2000/07/27 03:36:08  radphi
  Changes made by gevmage@urs3
  added proper updating of constants from map

  Revision 1.3  2000/06/25 08:05:47  radphi
  Changes made by radphi@urs3
  adding code to manipulate maps; not quite finished.

  Revision 1.2  2000/06/18 06:38:50  radphi
  Changes made by radphi@ifarms2
  now contains code to get BSD ccs

  Revision 1.1  1997/06/02 20:45:48  radphi
  Initial revision by lfcrob@dustbunny
  Basic cc management


  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calibration.h>
#include <map_manager.h>
#include <itypes.h>

static const char rcsid[]="$Id: calibration.c,v 1.7 2001/11/10 16:09:03 radphi Exp $";

static char *mapFile=NULL;
static char *map = "/maps/calibration.map";

int ccParameterLoad(int parameter_name, int runNo,
		    float *parameter, int *actualRun){
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"ccLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  switch(parameter_name){
  case PARAMETER_EPSILON:
    sprintf(str,"epsilon");        
    if(map_get_float(mapFile,"LGD",str,1,parameter,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load epsilon for run %d\n",runNo);
      return(1);
    }
    break;
  case PARAMETER_BETA:
    sprintf(str,"beta");        
    if(map_get_float(mapFile,"LGD",str,1,parameter,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load beta for run %d\n",runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"Unknown calibration parameter: %#06.4x\n",
	    parameter_name);
    return(1);
    break;
  }
  return(0);
}

int lgdTuneLoad(int runNo, int nChannels, float *cc, int *actualRun){
  char *dir; 
  char str[20];
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"ccLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  sprintf(str,"tune_%d",nChannels);        
  if(map_get_float(mapFile,"LGD",str,nChannels,cc,runNo,actualRun)){
    fprintf(stderr,"lgdTuneLoad: Failed to load tuned lgd constants for run %d\n",runNo);
    return(1);
  }
  return(0);
}

int ccLoad(int detector, int runNo, int nChannels, float *cc, int *actualRun)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"ccLoad: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  switch(detector){
  case CC_LGD:
    sprintf(str,"cc_%d",nChannels);        
    if(map_get_float(mapFile,"LGD",str,nChannels,cc,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load LGD cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_CPV:
    sprintf(str,"cc_%d",nChannels);        
    if(map_get_float(mapFile,"CPV",str,nChannels,cc,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load CPV cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_UPV:
    sprintf(str,"cc_%d",nChannels);        
    if(map_get_float(mapFile,"UPV",str,nChannels,cc,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load UPV cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_RPD:
    sprintf(str,"cc_%d",nChannels);        
    if(map_get_float(mapFile,"RPD",str,nChannels,cc,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load RPD cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_ESUM:
    sprintf(str,"cc_%d",nChannels);        
    if(map_get_float(mapFile,"ESUM",str,nChannels,cc,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load ESUM cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_BSD:
    sprintf(str,"cc_%d",nChannels);
    if(map_get_float(mapFile,"BSD",str,nChannels,cc,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load BSD cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_BGV:
    sprintf(str,"cc_%d",nChannels);
    if(map_get_float(mapFile,"BGV",str,nChannels,cc,runNo,actualRun)){
      fprintf(stderr,"ccLoad: Failed to load BGV cc's for run %d\n",runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"ccLoad : Unknown detector (0x%x)\n",detector);
    return(1);
  }  
  return(0);
}


int ccSave(int detector, int runNo, int nChannels, float *cc)
{
  char *dir; 
  char str[20];
  
  if(!mapFile){
    if((dir = getenv("RADPHI_CONFIG")) == NULL){
      fprintf(stderr,"ccSave: $RADPHI_CONFIG must be set\n");
      return(1);
    }
    mapFile = malloc(strlen(dir)+strlen(map)+1);
    sprintf(mapFile,"%s%s",dir,map);
  }
  runNo = abs(runNo);
  switch(detector){
  case CC_LGD:
    sprintf(str,"cc_%d",nChannels);        
    if(map_put_float(mapFile,"LGD",str,nChannels,cc,runNo)){
      fprintf(stderr,"ccSave: Failed to save LGD cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_CPV:
    sprintf(str,"cc_%d",nChannels);        
    if(map_put_float(mapFile,"CPV",str,nChannels,cc,runNo)){
      fprintf(stderr,"ccSave: Failed to save CPV cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_UPV:
    sprintf(str,"cc_%d",nChannels);        
    if(map_put_float(mapFile,"UPV",str,nChannels,cc,runNo)){
      fprintf(stderr,"ccSave: Failed to save UPV cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_RPD:
    sprintf(str,"cc_%d",nChannels);        
    if(map_put_float(mapFile,"RPD",str,nChannels,cc,runNo)){
      fprintf(stderr,"ccSave: Failed to save RPD cc's for run %d\n",runNo);
      return(1);
    }
    break;
  case CC_ESUM: 
    sprintf(str,"cc_%d",nChannels);        
    if(map_put_float(mapFile,"ESUM",str,nChannels,cc,runNo)){
      fprintf(stderr,"ccSave: Failed to save ESUM cc's for run %d\n",runNo);
      return(1);
    }
    break;
  default:
    fprintf(stderr,"ccSave : Unknown detector (0x%x)\n",detector);
    return(1);
  }  
  return(0);
}

