/*
  $Log: cern_minuit.c,v $
  Revision 1.2  2004/05/06 20:04:34  radphi

  D.K.

  Revision 1.1  1997/05/03 20:37:16  radphi
  Initial revision by lfcrob@dustbunny
  Minuit wrapper routines

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cernlib.h>

static const char rcsid[] = "$Id: cern_minuit.c,v 1.2 2004/05/06 20:04:34 radphi Exp $";

void mninit_(int *,int *,int *);
void mnparm_(int *,char *,double *,double *,double *,double*,int *,unsigned long);
void mnseti_(char *,unsigned long);
void mnexcm_(void (*fcn)(),char *,double *,int *,int *,void (*futil)(),unsigned long);
void mnpout_(int *,char *,double *,double *,double *,double *,int *,unsigned long);
void mnerrs_(int *,double *,double *,double *,double *);
void mnemat_(double **,int *);
void mnstat_(double *,double *, double *,int *, int *, int *);
void mninit(int i, int j, int k)

{
  mninit_(&i,&j,&k);
}

void mnparm(int index,char *name,double initialValue,double initialStep,
	    double min,double max,int flag)
{
  mnparm_(&index,name,&initialValue,&initialStep,&min,&max,&flag,strlen(name));
}

void mnseti(char *name)
{
  mnseti_(name,strlen(name));
}

void mnexcm(void (*fcn)(),char *command,double *arglist,int narg,int erflag,void (*futil)())
     
{
  mnexcm_(fcn,command,arglist,&narg,&erflag,futil,strlen(command));
}

void mnpout(int index,char *name,double *value, double *error,double *min,double *max,int *internalIndex)

{
  int length;

  mnpout_(&index,name,value,error,min,max,internalIndex,80);
}

void mnerrs(int index,double *plus,double *minus,double *para,double *globcc)
{
  mnerrs_(&index,plus,minus,para,globcc);
}

void mnemat(double **covar,int ndim)

{
  mnemat_(covar,&ndim);
}
void mnstat(double *fmin,double *fedm,double *errdef,int *npari,int *nparx,int *istat){
	  mnstat_(fmin,fedm,errdef,npari,nparx,istat);
}
