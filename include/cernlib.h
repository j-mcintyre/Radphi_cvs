
#ifndef CERNLIB_H_INCLUDED
#define CERNLIB_H_INCLUDED

#endif

/*
  IMPORTANT!!!!!

  You must declare the paw common block:

  struct hbk{float ia[800000];} pawc_;
*/

void save_histo(int, char *);



void hlimit(int size);
void hlimap(int size, char *name);
void hrput(int id, char *filename, char*command);
void hrget(int id, char *filename, char*command);
void hidopt(int id, char *command);



void hbookn(int,char *,int,char *,int,char **);
void hbook2(int,char *,int,float,float,int,float,float,float);
void hbook1(int,char*,int,float,float,float);
void hfill(int,float,float,float);
void hf2(int,float,float,float);
void hf1(int,float,float);
void hfn(int,float *);
void hunpak(int,float *,char *,int);
void hpak(int,float *);
void hldir(char *,char *);
void hmdir(char *,char *);
void hcdir(char *,char *);
void hropen(int,char *,char *,char *,int *,int *);
void hrout(int,int,char *);
void hrin(int, int, int);
void hrend(char *);
void hreset(int,char *);
void hbnt(int id,char*chtitle,char*chopt);
void hbname(int id,char*chblok,void*variable,char*chform);
void hfnt(int id);
void hfntb(int id,char* chblok);
void hprnt(int id);
void hfitga(int id, float *constant, float *average, float *sigma,float *chisq,float errors[3]);
void hplzon(int,int,int,char *);
void hplcon(int,int,int);
void hplego(int,float,float);
void iuwk(int,int);
void hplot(int,char *,char *,int);
void hplint(int);
void hgiven(int id, char *chtitl, int *nvar, char *chtag[], float rlow[], float rhigh[]);
void hgn(int id, int idn, int idnevt, float x[], int *ierror);
void hgnpar(int id, char *chrout);
void hgnf(int id, int idnevt, float x[], int *ierror);
void hnoent(int id, int *noent);
void hgnt(int id, int idnevt, int *ierr);
void hgntb(int id, char *chblok, int irow, int *ierr);
void hgntv(int id, char *chvar[], int nvar, int irow, int *ierr);
void hgntf(int id, int irow, int *ierr);
void hprntu(int id);

void fitGauss(int id, int lower, int upper, float par[3], float *chisq,
              float errors[3]);


/* For minuit */

void mninit(int i, int j, int k);
void mnparm(int index,char *name,double initialValue,double initialStep,
            double min,double max,int flag);
void mnseti(char *);
void mnexcm(void (*fcn)(),char *command,double *arglist,int narg,int erflag,void (*futil)());
void mnpout(int index,char *name,double *value, double *error,double *min,double *max,int *internalIndex);
void mnerrs(int index,double *plus,double *minus,double *para,double *globcc);
void mnemat(double **covar,int ndim);
void mnstat(double *fmin,double *fedm,double *errdef,int *npari,int *nparx,int *istat);


