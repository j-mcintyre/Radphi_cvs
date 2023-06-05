#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ReadUserIntensity.h>

int ReadIntensityFile(char *filename,float **I, int *nMassBins, float *massLow, float *massBinSize,
		      int *nCthBins, int *nPhiBins, float *Scale)
{
  FILE *fptr=NULL;
  int fd;
  int ret=0;
  int dummy;
  if(!(fptr=fopen(filename,"r"))){
    fprintf(stderr,"\nReadIntensityFile: unable to open file %s ... returning failure",filename);
    ret=1;
  }
  else{
    fprintf(stderr,"\nReadIntensityFile: opened file %s ...",filename);
    fd=fileno(fptr);
    
    fprintf(stderr,"\n=================== \nreading intensity array:");    
    dummy=read(fd,nMassBins,sizeof(int));
    fprintf(stderr,"\n %d mass bins",*nMassBins);
    dummy=read(fd,massLow,sizeof(float));
    fprintf(stderr,"\n %f low mass",*massLow);
    dummy=read(fd,massBinSize,sizeof(float));
    fprintf(stderr,"\n %f mass binsize",*massBinSize);
    dummy=read(fd,nCthBins,sizeof(int));
    fprintf(stderr,"\n %d bins in costheta",*nCthBins);
    dummy=read(fd,nPhiBins,sizeof(int));
    fprintf(stderr,"\n %d bins in phi",*nPhiBins);
    if(Scale!=NULL){
      Scale =(float*)malloc(*nMassBins*sizeof(float));
    }
    else{
      fprintf(stderr,"ReadIntensityFile: non-NULL pointer Scale detected, will free, re-malloc, and continue...cross fingers.");
      free(Scale);
      Scale=NULL;
      Scale =(float*)malloc(*nMassBins*sizeof(float));
    }
    dummy=read(fd,Scale,*nMassBins*sizeof(float));
    fprintf(stderr,"\n scale factor array read");

    if(*I==NULL){
      /* then it's ok to malloc, else take exception */
      *I=(float*)malloc((*nMassBins)*(*nCthBins)*(*nPhiBins)*sizeof(float));
    }
    else{
      fprintf(stderr,"ReadIntensityFile: non-NULL pointer I detected, will free, re-malloc, and continue running \nPLEASE FIX THIS");
      free(I);
      *I=(float*)malloc((*nMassBins)*(*nCthBins)*(*nPhiBins)*sizeof(float));
    }
    dummy=read(fd,*I,(*nMassBins)*(*nCthBins)*(*nPhiBins)*sizeof(float));
  }
  fprintf(stderr,"\n===== done reading ============\n");    
  return (ret);
}
  
  
