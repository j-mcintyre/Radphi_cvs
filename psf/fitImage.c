#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <dataIO.h>
#include <disData.h>

#include <fileIO.h>
#include <eventIO.h>
#include <psf.h>
#include <cernlib.h>

/* private data */

#define SQR(x) ((x)*(x))
#define LGD_MATRIX_ROWS 28
#define LGD_MATRIX_COLS 28
#define FIT_ROW_OFFSET 3
#define FIT_COL_OFFSET 3
#define LGD_PLANE_DISTANCE 120
#define MAX_ENERGY_SLICES 100 
#define SIGNIFICANT_ROW_OFFSET 3
#define SIGNIFICANT_COL_OFFSET 3
#define MOMENT_ORDER 4

/* public function */

double fitgetchi2_(int *npar, double *pars, int *ndof);
int fithf_(int *npar, double *pars);
int fittest_(void);
int fitgetmom_(void);

void fitgetpar0_(int *npar,double *pars);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Dummy call to load the module
*/
void fitImage(void) {}           /* c entry point */
void fitimage_(void) {}          /* fortran entry */

/* private function to find out moments, central moments 
* of x and y coordinates of hits within a shower image. 
*/ 
static int fitGetMomentXY(int nMom, float *xMom, float *yMom);
static int fitGetCentralMoment(int nMom, float *aMom, float *cMom);

/* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User call-able method to get chi2 for the given
* set of points and fitting function
*/
double fitgetchi2_(int *npar, double *pars, int *ndof) {
 
 int ret;
 int ch,row,col;
 int rowMax,colMax;
 float Emax;
 double D,E,Sigma2;
 double chi2=0;
 int ndf=0;

 if ( (ret=psfSetPars(*npar,pars)) != 0) {
   return -1.;
 } 
 
 eventGetEmax(&ch);
 eventGetLgdRC(ch,&rowMax,&colMax); 
 

/* Change the skim to avoid "zero" bug: fit around max energy
   but do not drop zero-data from  fit! In addition, we will always 
   have some error values thus no need to vory about Sigma2=0! 
*/
 for(row=(rowMax-FIT_ROW_OFFSET);row<=(rowMax+FIT_ROW_OFFSET);row++) {
   for(col=(colMax-FIT_COL_OFFSET);col<=(colMax+FIT_COL_OFFSET);col++) {
     if ( row>=0 && row<LGD_MATRIX_ROWS &&
          col>=0 && col<LGD_MATRIX_COLS ) {
       ch = row*LGD_MATRIX_ROWS + col;
       
       D=eventGetErc(row,col);
       E=psfGetE(row,col);
       Sigma2=psfGetSigma2(row,col);
       chi2 += SQR(D-E)/Sigma2;
       ndf++;
     }
/* Fill histograms for testing (do not during fitting)!!! 
     hf2(200,(float) col,(float) row,D);
     hf2(201,(float) col,(float) row,E);
     hf2(202,(float) col,(float) row,Sigma2); 
*/
   }
 }
 *ndof = ndf;
 return chi2;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User call-able method to get initial parameter values from data
*/
void fitgetpar0_(int *npar,double *pars) {
  
  float Emax;
  int ch,row,col;
  float X0,Y0;

  Emax = eventGetEmax(&ch);
  eventGetLgdRC(ch,&row,&col);
  eventGetLgdXY(row,col, &X0, &Y0);
  pars[0] = (double) eventGetEtot();
  pars[1] = (double) X0;
  pars[2] = (double) Y0;
  return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User call-able method to test chi2 func 
*/
int fittest_(void){
  
  int ret;
  int index;

  float Emax;
  int ch,row,col;
  float X0,Y0;

  int nparam=4;
  int nDoF;
  double pars[nparam];
  double chi2;

  if ((ret=eventGetNext()) < 0) {
    fprintf(stderr,"ERROR - cannot get an event!");
    return -1;
  }
  fprintf(stderr," Processing %d event, E %f, ...\n",ret,eventGetEtot());

  hf1(20,eventGetEtot(),1.) ;
  fitgetpar0_(&nparam,pars);
  chi2 = fitgetchi2_(&nparam,pars,&nDoF);
  fprintf(stderr,"Chi2: %f, Ndof: %d, (Chi2/Ndof): %f\n",chi2,nDoF,chi2/nDoF);
  return 0;  
}

/* * * * * * * * * * * * * * * * * * * * * * * *
* User call-able function for histogram filling.
* Fortran interface!
*/ 
int fithf_(int *npar, double *pars) {

  int ret;
  int i;
  int ch,row,col;
  int chCent=0, rowCent=0, colCent=0;
  float Etot,Emax;
  float Ephot,Emc,Eref;
  double D,E,Sigma,R;
  double chi2=0;
  float Z0 = LGD_PLANE_DISTANCE;
  float EI,E0,E2;
  float Efit,X,Y,theFit,phiFit;
  float px,py,pz,theta,phi;
  float pxMC,pyMC,pzMC,theMC,phiMC;
  float sinx,siny;
  float sinxf,sinyf;
  
  float Ebin;
  int nbins = 9;
  static float Ebond[] = 
      {0., 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 2.8, 3.6, 6.0};
  
  int id1,id2,id3,id8,id9;
  int id11,id81;
  int rNgb,cNgb;

  if ( (ret=psfSetPars(*npar,pars)) != 0) {
   return -1.;
  }
  Etot=eventGetLgdEtot();
  Emax=eventGetLgdEmax(&ch);
  hf1(20,Etot,1.);
  hf1(21,Emax,1.);

  Etot=eventGetEtot();
  Emax=eventGetEmax(&ch);
  Ephot=eventphoton_(&px,&py,&pz);
  theta = atan2(sqrt(SQR(px)+SQR(py)),pz);
  phi = atan2(py,px);

  hf1(22,Etot,1.);
  hf1(23,Emax,1.);
  hf1(24,Ephot,1.);
//  hf2(25,Ephot,Emax,1.);
  hf1(26,theta,1.);
  hf1(27,phi,1.);
  
  Efit = pars[0];
  X = pars[1];
  Y = pars[2];
  R = pars[3];
  theFit = atan2(sqrt(X*X+Y*Y),Z0);
  phiFit = atan2(Y,X);
  sinxf = sin(theFit)*cos(phiFit);
  sinyf = sin(theFit)*sin(phiFit);

  hf1(300,Efit,1.);
  hf1(301,X,1.);
  hf1(302,Y,1.);
  hf1(303,R,1.);
  hf1(306,theFit,1.);
  hf1(307,phiFit,1.);
  
  hfill(400,Ephot,Efit,1.);
  hfill(401,Ephot,R,1.);
  hfill(402,Ephot,theFit,1.);
  hfill(403,Ephot,phiFit,1.);
  hfill(500,Ephot,Efit,1.);
  hfill(501,Ephot,R,1.);
  hfill(502,Ephot,theFit,1.);
  hfill(503,Ephot,phiFit,1.);

  if (eventmc_() > 0) {
    Emc = mcphoton_(&pxMC,&pyMC,&pzMC);
    Eref = Emc;
    theMC = atan2(sqrt(SQR(pxMC)+SQR(pyMC)),pzMC);
    phiMC = atan2(pyMC,pxMC);
    hf1(25,Emc,1.);
    hf1(28,theMC,1.);
    hf1(29,phiMC,1.);
    hfill(404,Emc,Efit,1.);
    hfill(405,Emc,R,1.);
    hfill(406,Emc,theFit,1.);
    hfill(407,Emc,phiFit,1.);
    hfill(408,Emc,Etot,1.);
    hfill(504,Emc,Efit,1.);
    hfill(505,Emc,R,1.);
    hfill(506,Emc,theFit,1.);
    hfill(507,Emc,phiFit,1.);
    hfill(508,Emc,Etot,1.);
    hfill(509,Emc,Ephot,1.);

    sinx = pxMC/Emc;
    siny = pyMC/Emc;
    if ( (fabs(sinx-sinxf)>0.1) || (fabs(siny-sinyf)>0.1 ) ) {
      fprintf(stderr,"MC: E %f\t theta %f\t phi %f\n",
                      Emc,theMC,phiMC); 
      fprintf(stderr,"Fit: E %f\t theta %f\t phi %f\n",
                      Efit,theFit,phiFit);  
    } 
  }
  else {
    Eref = Ephot;
    sinx = px/Ephot;
    siny = py/Ephot;
    if ( (fabs(sinx-sinxf)>0.1) || (fabs(siny-sinyf)>0.1 ) ) {
      fprintf(stderr,"Phot: E %f\t theta %f\t phi %f\n",
                      Ephot,theta,phi); 
      fprintf(stderr,"Fit: E %f\t theta %f\t phi %f\n",
                      Efit,theFit,phiFit);  
    } 
  }
/* Histogram ID's for data slced in energy 
* Convenction: id1 - data, id2 - fit, id3 - errors, 
*              id8,id9, - data,fit bin by bin.   
* For later normalization the number of hits is also recorded using
* idX1 for hits in slice X.*/  
  id1=1000;
  id2=2000;
  id3=3000;
  id8=8000;
  id9=9000;

  id11=1100;
  id81=8100;
/* fill LGD hits, PSF hits, errors, and yields of neighbouring blocks 
* around central image channel for all energies!
*/   
  Emax=eventGetEmax(&ch);
  eventGetLgdRC(ch,&rowCent,&colCent);
  chCent=ch;
  rNgb = SIGNIFICANT_ROW_OFFSET;
  cNgb = SIGNIFICANT_COL_OFFSET;

 for(row=0;row<LGD_MATRIX_ROWS;row++){
   for(col=0;col<LGD_MATRIX_COLS;col++){
     D=eventGetErc(row,col);
     E=psfGetE(row,col);
     Sigma=sqrt(psfGetSigma2(row,col));
     hf2(200,(float) col,(float) row,D);
     hf2(201,(float) col,(float) row,E);
     hf2(202,(float) col,(float) row,Sigma);
   }	  
 } 
 ch=0;
 for(row=(rowCent-rNgb);row<=(rowCent+rNgb);row++) {
   for(col=(colCent-cNgb);col<=(colCent+cNgb);col++) {
     if ( row>=0 && row<LGD_MATRIX_ROWS &&
          col>=0 && col<LGD_MATRIX_COLS ) {
       D=eventGetErc(row,col);
       E=psfGetE(row,col); 
       hf1(208,(float) ch,D);
       hf1(209,(float) ch,E);
       ch++;
     }  
   }
 }
/* fill LGD hits, PSF hits, errors, and yields of neighbouring blocks 
* around maxChannel for all energy slices!
*/   
 for(i=0;i<(nbins-1);i++){
   Ebin=-1;	 
   if ((Ebond[i]<Eref)&&(Eref<Ebond[i+1])) {
     Ebin = 0.5*(Ebond[i]+Ebond[i+1]);	   
     hf1(33,Ebin,1.);
     for(row=0;row<LGD_MATRIX_ROWS;row++){
       for(col=0;col<LGD_MATRIX_COLS;col++){
         D=eventGetErc(row,col);
         E=psfGetE(row,col);
         Sigma=sqrt(psfGetSigma2(row,col));
         hf2(id1+i,(float) col,(float) row,D);
         hf2(id2+i,(float) col,(float) row,E);
         hf2(id3+i,(float) col,(float) row,Sigma);	   
         hf2(id11+i,(float) col,(float) row,1.);
       }
     }
     ch=0;
     for(row=(rowCent-rNgb);row<=(rowCent+rNgb);row++) {
       for(col=(colCent-cNgb);col<=(colCent+cNgb);col++) {
         if ( row>=0 && row<LGD_MATRIX_ROWS &&
              col>=0 && col<LGD_MATRIX_COLS ) {
           ch++;
	   D=eventGetErc(row,col);
	   E=psfGetE(row,col); 
	   hf1(id8+i,(float) ch,D);
	   hf1(id9+i,(float) ch,E);
	   hf1(id81+i,(float) ch,1.);
	 }  
       }
     }
   }
 }
  return 0;
}


/* User call-able mothod to plot moments, central moments as well as
*  Variance, Skewness, and Kurtosis of X and Y hit coordinates within the
* shower image.
*/
int fitgetmom_(void) {

  int nMom=MOMENT_ORDER+1;
  int ret;
  float aMomX[nMom];
  float aMomY[nMom];
  float cMomX[nMom];
  float cMomY[nMom];
  float VarianceX=0.;
  float SkewnessX=0.;
  float KurtosisX=0.;
  float VarianceY=0.;
  float SkewnessY=0.;
  float KurtosisY=0.;
  float eps=1.e-6;
  float M0=0;

  memset(aMomX,0.,sizeof(float)*nMom);
  memset(aMomY,0.,sizeof(float)*nMom);
  memset(cMomX,0.,sizeof(float)*nMom);
  memset(cMomY,0.,sizeof(float)*nMom);
  
  ret=fitGetMomentXY(nMom,aMomX,aMomY);
  if(ret<0) return -1;
  
  M0 = aMomX[0];
  ret=fitGetCentralMoment(nMom,aMomX,cMomX);
  if(ret<0) return -1;
  ret=fitGetCentralMoment(nMom,aMomY,cMomY);
  if(ret<0) return -1;
	  
  VarianceX = cMomX[2];
  SkewnessX = cMomX[3]/(pow(VarianceX,1.5)+eps);
  KurtosisX = cMomX[4]/(pow(VarianceX,2)+eps) - 3. ;
  VarianceY = cMomY[2];
  SkewnessY = cMomY[3]/(pow(VarianceY,1.5)+eps);
  KurtosisY = cMomY[4]/(pow(VarianceY,2)+eps) - 3. ;

  if (M0>1.0) { 
    hfill(350,aMomX[1]/M0,0.,1.);
    hfill(351,sqrt(VarianceX),0.,1.);
    hfill(352,SkewnessX,0.,1.);
    hfill(353,KurtosisX,0.,1.);
  
    hfill(360,aMomY[1]/M0,0.,1.);
    hfill(361,sqrt(VarianceY),0.,1.);
    hfill(362,SkewnessY,0.,1.);
    hfill(363,KurtosisY,0.,1.);
  }  
  return 0;  
}


static int fitGetMomentXY(int nMom, float *xMom, float *yMom) {

  int index;
  int ch,row,col;
  int rowMax,colMax;
  float E;
  float X,Y;

  E=eventGetEmax(&ch);
  eventGetLgdRC(ch,&rowMax,&colMax);
  if (E==0) return -1;

  for(row=(rowMax-FIT_ROW_OFFSET);row<=(rowMax+FIT_ROW_OFFSET);row++) {
    for(col=(colMax-FIT_COL_OFFSET);col<=(colMax+FIT_COL_OFFSET);col++) {
      if ( row>=0 && row<LGD_MATRIX_ROWS &&
           col>=0 && col<LGD_MATRIX_COLS ) {
        E=eventGetErc(row,col);
        eventGetLgdXY(row,col,&X,&Y);
        for (index=0;index<nMom;index++) {
          xMom[index] += E*pow(X,index);
          yMom[index] += E*pow(Y,index);
        }
      }
    }
  }

}

/* Private method to calculate Central Moments (cMom) from the
* given Algebraic Moments (aMom) of order nMom.
*/  
static int fitGetCentralMoment(int nMom, float *aMom, float *cMom) {

  float eps=1.e-6;
  int i;
  float cm[nMom];
  
  memset(cm,0,sizeof(float)*nMom);
  
  switch (nMom-1){
    case 4:
      cm[4] = (aMom[4]*pow(aMom[0],3)-4.*aMom[3]*aMom[1]*pow(aMom[0],2)
		                     +6.*aMom[2]*pow(aMom[1],2)*aMom[0]
	                             -3.*pow(aMom[1],4))
	                                    /(pow(aMom[0],4)+eps);
    case 3:
      cm[3] = (aMom[3]*pow(aMom[0],2)-3.*aMom[2]*aMom[1]*aMom[0]+2.*pow(aMom[1],3))
                                            /(pow(aMom[0],3)+eps);
    case 2:
      cm[2] = (aMom[2]*aMom[0]-pow(aMom[1],2))/(pow(aMom[0],2)+eps) ;
    case 1:
      cm[1] = 0.;
    case 0:
      cm[0] = 1.;
      break;
    default:
       fprintf(stderr,"Warning - central moment of order greater than %d",nMom-1);
       fprintf(stderr,"          not defined (yet)\n");
       return -1;
  }

  for(i=0; i<nMom; i++) {
    cMom[i]=cm[i];
  }
  return 0; 
}

