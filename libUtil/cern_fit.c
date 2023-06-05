#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cernlib.h>

/*
 $Log: cern_fit.c,v $
 Revision 1.1  1996/11/21 12:39:18  lfcrob
 Initial revision

*/

static const char rcsid[]="$Id: cern_fit.c,v 1.1 1996/11/21 12:39:18 lfcrob Exp $";

void hfithn_(int *,char *,char *,int *,float *,float *,float *,float *,float *,float *,unsigned long,unsigned long);




struct fit{int iquest[100];} quest_;


void fitGauss(int id, int lower, int upper,float par[3],float *chisq,float errors[3])
{

  char *chfun = "G";
  char *chopt = "FQR";
  int np=3;
  float step[3];
  float pmin[3];
  float pmax[3];
  


  quest_.iquest[10] = lower;
  quest_.iquest[11] = upper;

  hfithn_(&id,chfun,chopt,&np,par,step,pmin,pmax,errors,chisq,strlen(chfun),strlen(chopt));
}
