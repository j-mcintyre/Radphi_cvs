#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calibration.h>
#include <map_manager.h>
#include <detector_params.h>
#include <itypes.h>

#define LENGTH_OF_ARRAY 38
#define LENGTH_OF_PATHNAME 256

enum {SUCCESS, FAIL};

int main(int nargs, char *argv[]){
  int i;
  FILE *fptr;
  char *filename; 
  char *radphi_config;
  char response;
  int reval = SUCCESS;
  int runnumber; 
  int tagnum[LENGTH_OF_ARRAY];
  float mean[LENGTH_OF_ARRAY];
  float sigma[LENGTH_OF_ARRAY];
  runnumber = 0;
  
  if ((nargs != 3)){
    fprintf(stderr,"need two command line arguments,\n");  
    fprintf(stderr,"A run number and filename.\n");
    fprintf(stderr,"Example: data_put 8012 tagger_mean8012.dat\n"); 
    exit (1);
  }
  runnumber = atoi(argv[1]);
  filename = argv[2];
  printf("Run Number = %d\n", runnumber);
  
  if ((fptr = fopen(filename, "r")) == NULL){
    printf("Cannot open %s.\n", filename);
    reval = FAIL;
    exit(1);
  } 
  else { 
    /* file is open, read tagger number, mean time and time width*/
    
    int i;
    
    /* loop over number of tagger channels (18 left, 18 right)*/
    
    for (i=0; i<LENGTH_OF_ARRAY; i++){ 
      fscanf(fptr, "%d    %f    %f\n", &tagnum[i], &mean[i], &sigma[i]);
      /*      printf("%d\t%.2f\t%.2f\n",tagnum[i], mean[i], sigma[i]);*/
      
    }
    fclose(fptr);
  }
  
  /* Save data into map tdc_offset_tag*/
  
  fprintf(stderr,"I have finished parsing your data file.\n");
  fprintf(stderr,"The data that is going to be sent to the map is:\n");
  for(i=0;i<LENGTH_OF_ARRAY;i++){
    fprintf(stderr,"index: %3d value: %f\n",i,mean[i]);
  }
  
  fprintf(stderr,"Do you want to save this data to the map? ");
  response=getchar();
  if(!(response=='y'||response=='Y')){
    fprintf(stderr,"Aborting Upload.\n");
    exit(1);
  }
  radphi_config=getenv("RADPHI_CONFIG");
  if(!strcmp(radphi_config,"/home/radphi/config")){
    fprintf(stderr,"\n\nBZZZT!   BZZZT!   BZZZT!\n");
    fprintf(stderr,"WARNING!  WARNING!  You have invoked data_put to\n");
    fprintf(stderr,"upload constants into the map while your\n");
    fprintf(stderr,"RADPHI_CONFIG was pointing at the global radphi\n");
    fprintf(stderr,"config area.\n\n");
    fprintf(stderr,"Please change this and try again.\n");
    fprintf(stderr,"An example would be (enter at the command line):\n");
    fprintf(stderr,"setenv RADPHI_CONFIG ~your_username/config\n");
    exit(1);
  }
  
  if(tdcOffsetSave(TDC_OFFSET_TAG,runnumber,LENGTH_OF_ARRAY,mean)){
    fprintf(stderr,"error saving sample data into TAG tdc offsetmap.\n");
    exit(FAIL);
  }
  
  fprintf(stderr,"Save successful.\n");
  return reval;
}       




