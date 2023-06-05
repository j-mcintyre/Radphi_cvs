#include <stdio.h>
#include <unistd.h>

#ifndef _READ_USER_INTENSITY_INCLUDED_
#define _READ_USER_INTENSITY_INCLUDED_

int  ReadIntensityFile(char *filename,float **I,int *nMassBins, 
		       float *massLow, float *massBinSize,	
		       int *nCthBins, int *nPhiBins, float *Scale);


#endif
