#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calibration.h>
#include <map_manager.h>
#include <detector_params.h>
#include <itypes.h>
#include <makeHits.h>

int main(int nargs, char args[]){
  int i;
  float sample_array[N_TAG_TDC_OFFSET];
  int run_number, taggn, mean;
  FILE *fptr;
  
  
  for(i=0;i<N_TAG_TDC_OFFSET;i++){
    sample_array[i]=.1+(((float)i)*(0.3));
     
    enum {SUCCESS, FAIL};
    void DataWrite (File *fout);
    void DataRead(File *Tagger_mean[run_number].dat);
   
  {
    char []= "strnum.mix";
    int reval = SUCCESS;
        

    if ((fptr = fopen(, "W+")) == NULL){
       reval = ErrorMsg(filenam);
    }  else {
        
    
        fscanf(FILE *fptr, "%t %d %t %d", taggn, mean);   
    
   
    return reval;
  }
  
  /* function definition*/
void DataWrite (File *TDC_OFFSET_TAG)
  {
    int i;
    char  means[MAX_NUM] = {
     


 /*}
   if(tdcOffsetSave(TDC_OFFSET_TAG,run_number,N_TAG_TDC_OFFSET,sample_array)){
    fprintf(stderr,"error saving sample data into TAG tdc offsetmap.\n");
    exit(1);
  }
  fprintf(stderr,"Save successful.\n");

  }*/



