#include <stdio.h>
#include <stdlib.h>
#include <makeHits.h>
#include <map_editing.h>
/* #define DEBUGGING 1 */

int setup_global_utility_structures(void){
  int i;
#ifdef DEBUGGING 
  fprintf(stderr,"begin setup_global_utility_structures.\n");
#endif
  n_channels_in_detector[GLOBAL_DETECTOR_NUMBER_TAG]=
    N_TAG_CHANNELS;
  n_channels_in_detector[GLOBAL_DETECTOR_NUMBER_UPV]=
    N_UPV_CHANNELS;
  n_channels_in_detector[GLOBAL_DETECTOR_NUMBER_BSD]=
    N_BSD_CHANNELS;
  n_channels_in_detector[GLOBAL_DETECTOR_NUMBER_BGV]=
    N_BGV_CHANNELS;
  n_channels_in_detector[GLOBAL_DETECTOR_NUMBER_CPV]=
    N_CPV_CHANNELS;
  n_channels_in_detector[GLOBAL_DETECTOR_NUMBER_LGD]=
    lgd_nChannels;

#ifdef DEBUGGING 
  fprintf(stderr,"n_channels_in_detector_array assigned.\n");
#endif

  for(i=0;i<N_GLOBAL_DETECTOR_NUMBERS;i++){
    if(detector_name[i]==NULL){
      detector_name[i]=
	(char*)malloc(sizeof(char)*DETECTOR_NAME_STRING_LENGTH);
      sprintf(detector_name[i],"BAD");
    }
  }
  if(N_GLOBAL_DETECTOR_NUMBERS<6){
    fprintf(stderr,"setup_global_utility_structures:\n");
    fprintf(stderr,"N_GLOBAL_DETECTOR_NUMBERS<6!!  ARGH!\n");
    return(1);
  }
  sprintf(detector_name[GLOBAL_DETECTOR_NUMBER_TAG],"TAG");
  sprintf(detector_name[GLOBAL_DETECTOR_NUMBER_UPV],"UPV");
  sprintf(detector_name[GLOBAL_DETECTOR_NUMBER_BSD],"BSD");
  sprintf(detector_name[GLOBAL_DETECTOR_NUMBER_BGV],"BGV");
  sprintf(detector_name[GLOBAL_DETECTOR_NUMBER_CPV],"CPV");
  sprintf(detector_name[GLOBAL_DETECTOR_NUMBER_LGD],"LGD");
  
#ifdef DEBUGGING 
  fprintf(stderr,"end of setup_global_utility_structures.\n");
#endif
  
  return(0);
}
