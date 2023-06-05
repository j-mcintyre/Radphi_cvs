#include <stdio.h>
#include <stdlib.h>
enum {SUCCESS, FAIL,
      MAX = 38,
      NUM_COLUM = 3};

main(void) 
{
   FILE *fptr;
   char filename[]= "tagdata";
   int reval = SUCCESS;
   
   if ((fptr = fopen(filename, "r")) == NULL){
      printf("Cannot open %s.\n", filename);
      reval = FAIL;
   } 
   else { 
        int i; 
	int tagnum[MAX];
        float mean[MAX];
        float sigma[MAX];
        /*float tdc_offset_array[MAX][NUM_COLUM];*/
	/* float *tdc_offset;*/
	/* tdc_offset = tdc_offset_array;*/
	/*tdc_offset=malloc(sizeof(float)*MAX);*/
       for (i=0; i<MAX; i++){ 
         fscanf(fptr, "%d    %f    %f\n", &tagnum[i], &mean[i], &sigma[i]);
         printf("%d\t%.2f\t%.2f\n",tagnum[i], mean[i], sigma[i]);
                
       }
       fclose(fptr);
   }

    return reval;
}          
   
 



