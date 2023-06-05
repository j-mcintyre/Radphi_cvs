#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calibration.h>
#include <map_manager.h>
#include <detector_params.h>
#include <itypes.h>

#define LENGTH_OF_ARRAY 38

enum {SUCCESS, FAIL};

int main(int nargs, char *args[]){
  int i;
   FILE *fptr;
   char filename[]= "tagger_mean"[runnumber]".dat";
   int reval = SUCCESS;
   int runnumber = 0; 
   runnumber = atoi(argv[1]);
   int tagnum[LENGTH_OF_ARRAY];
   float mean[LENGTH_OF_ARRAY];
   float sigma[LENGTH_OF_ARRAY];

/*mean=malloc(sizeof(float)*LENGTH_OF_ARRAY);*/

   /*open file "tagdata". Should be connected via a soft link to data file
     e.g. ln -s tagger_mean7871.dat tagdata */

   if ((fptr = fopen(filename, "r")) == NULL){
      printf("Cannot open %s.\n", filename);
      reval = FAIL;
   } 
   else { 
     /* file is open, read tagger number, mean time and time width*/

        int i;

	/* loop over number of tagger channels (18 left, 18 right)*/

        for (i=0; i<LENGTH_OF_ARRAY; i++){ 
         fscanf(fptr, "%d    %f    %f\n", &tagnum[i], &mean[i], &sigma[i]);
         printf("%d\t%.2f\t%.2f\n",tagnum[i], mean[i], sigma[i]);
                
       }
       fclose(fptr);
   }

   /* Save data into map tdc_offset_tag

  if(tdcOffsetSave(TDC_OFFSET_TAG,runnumber,LENGTH_OF_ARRAY,mean)){
    fprintf(stderr,"error saving sample data into TAG tdc offsetmap.\n");
    exit(FAIL);
    }*/
    
    fprintf(stderr,"Save successful.\n");
    return reval;
}       

