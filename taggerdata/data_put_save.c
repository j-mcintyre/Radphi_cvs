#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calibration.h>
#include <map_manager.h>
#include <detector_params.h>
#include <itypes.h>

#define LENGTH_OF_ARRAY 38

int main(int nargs, char args[]){
  int i;
  float *sample_array=NULL;
  
  sample_array=malloc(sizeof(float)*LENGTH_OF_ARRAY);
  
  for(i=0;i<LENGTH_OF_ARRAY;i++){
    sample_array[i]=.1+(((float)i)*(0.3));
  }
  if(tdcOffsetSave(TDC_OFFSET_TAG,5,LENGTH_OF_ARRAY,sample_array)){
    fprintf(stderr,"error saving sample data into TAG tdc offsetmap.\n");
    exit(1);
  }
  fprintf(stderr,"Save successful.\n");

}
