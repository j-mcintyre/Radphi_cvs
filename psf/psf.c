#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <eventIO.h>

#include <psf.h>

#define SQR(x) ((x)*(x))
#define MAX_PAR 50
#define BLOCK_SIZE 4

double *fpar=NULL;

/* public functions defined in psf.h */
/* c interface
double psfGetE(int row, int col);
double psfGetSigma2(int row, int col);
int psfSetPars(int npar, double *pars);

double psfGetEI(int row,int col);
double psf1Gss0(int row,int col);
double psf1Gss2(int row, int col);

 * fortran interface

double psffunc_(int row, int col);
double psfsigma2_(int row, int col);

double psf1Exp(int row, int col);
*/

/* private functions */

double gaussCylindric0(float deltaX, float deltaY, float E, float R);
double gaussCylindricI(float deltaX, float deltaY, float E, float R);
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* User call-able method to set internal parameter's space from
* the input vector
*/
int psfSetPars(int npar, double *pars) {

  int index;
  if (npar>MAX_PAR) {
    fprintf(stderr,"ERROR - too many pars.! \n");
    return -1;
  }
  if (fpar==NULL) {
    fpar = malloc(sizeof(double)*MAX_PAR);
    if(fpar==NULL) {
      fprintf(stderr,"ERROR - could not allocate parameter's buffer\n");
      return -1;
    }
  }

  for(index=0;index<npar;index++) {
    fpar[index] = pars[index];
/* Printout for debug 
    fprintf(stderr,"Par(%d): %f \n",index,*(fpar+index)); */
  }

  return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* User call-able method to obtain energy yield for the given coordinates
* from the look-up table. 
*/
double psfGetE(int row, int col){
  return psfGetEI(row,col); 
}

double psfGetEI(int row, int col){
  
   float x,y;
   int L = BLOCK_SIZE;
   double deltaX;
   double deltaY;
   double intX;
   double intY;
   double Etot = fpar[0];
   double X0 = fpar[1];
   double Y0 = fpar[2];
   double sigma1 = fabs(fpar[3]);

   eventGetLgdXY(row,col,&x,&y);

   deltaX = x-X0;
   deltaY=  y-Y0;

   intX = 0.5*(erf((deltaX+L/2)/(sqrt(2.)*sigma1))
               -erf((deltaX-L/2)/(sqrt(2.)*sigma1)));
 
   intY = 0.5*(erf((deltaY+L/2)/(sqrt(2.)*sigma1))
               -erf((deltaY-L/2)/(sqrt(2.)*sigma1)));
   return Etot*intX*intY;	/* returns bin integral of the Gaussian*/
}


double psf1Gss0(int row, int col){
   float x,y;
   int L2 = SQR(BLOCK_SIZE);
   double rho2,S2;
   double Norm1,psf0;
   double Etot = fpar[0];
   double X0 = fpar[1];
   double Y0 = fpar[2];
   double sigma1 = fpar[3];

   eventGetLgdXY(row,col,&x,&y);

   rho2 = SQR(x-X0) + SQR(y-Y0);
   S2 = SQR(sigma1);
   Norm1 = Etot/(2*3.14159*S2);
   psf0= L2*Norm1*exp(-0.5*rho2/S2);

   return psf0;  	/* returns 0-th integral expansion approx. */
}


double psf1Gss2(int row, int col) {

   float x,y;
   int L2 = SQR(BLOCK_SIZE);
   double rho2,S2;
   double Norm1,psf2;
   double correct2;
   double Etot = fpar[0];
   double X0 = fpar[1];
   double Y0 = fpar[2];
   double sigma1 = fpar[3];

   eventGetLgdXY(row,col,&x,&y);

   rho2 = SQR(x-X0) + SQR(y-Y0);
   S2 = SQR(sigma1);
   correct2 = (rho2-2.*S2)*L2/(24.*SQR(S2));
   Norm1 = Etot/(2*3.14159*S2);
   psf2= L2*Norm1*exp(-0.5*rho2/(S2))*(1.+correct2);
   
   return psf2;   /* returns bin 2-nd order integral expansion */
}

double psfGetSigma2(int row,int col) {

  double A=0.08;
  double B=5.01;
  double E=eventGetErc(row,col);
  double Etot=eventGetEtot();

  double err = A*sqrt(E)+B;
  return err*err;
}

double psf1Exp(int row, int col) {

   float x,y;
   int L2 = SQR(BLOCK_SIZE);
   double rho;
   double Norm1,psf1;
   double Etot = fpar[0];
   double X0 = fpar[1];
   double Y0 = fpar[2];
   double slope1 = fpar[3];

   eventGetLgdXY(row,col,&x,&y);

   rho = sqrt(SQR(x-X0) + SQR(y-Y0));
   Norm1 = Etot/(2*3.14159*SQR(slope1));
   psf1= L2*Norm1*exp(-rho/slope1);

   return psf1;   /* returns bin integral estimate */
}


double psffunc_(int row, int col){
  return psfGetE(row,col);
}

double psfsigma2_(int row, int col){
  return psfGetSigma2(row,col);
}

/* private functions */

/* private method to get the value of symmetric cylindrical Gaussian.
 * deltaX and deltaY define the radial distance from the center of
 * the Gaussian, E is overall normalization and R is radial width.
 */ 
double gaussCylindric0(float deltaX, float deltaY, float E, float R) {
 
   float rho2 = SQR(deltaX) + SQR(deltaY);
   float S2 = SQR(R);
   float N = E/(2.*M_PI*S2);

   return N*exp(-0.5*rho2/S2);

}

/* private method to get the value of symmetric cylindrical Gaussian
 * integrated over the area of one block. 
 */
double gaussCylindricI(float deltaX, float deltaY, float E, float R) {

   int L = BLOCK_SIZE;
   double intX,intY;

   intX = 0.5*(erf((deltaX+L/2)/(sqrt(2.)*R))
               -erf((deltaX-L/2)/(sqrt(2.)*R)));
 
   intY = 0.5*(erf((deltaY+L/2)/(sqrt(2.)*R))
               -erf((deltaY-L/2)/(sqrt(2.)*R)));
   return E*intX*intY;
}

/* private method to get the value of symmetric cylindrical Gaussian
 *
 */
float shiftCenter() {
}

double psf2Gss2(int row, int col) {

   float x,y;
   int L2 = SQR(BLOCK_SIZE);
   double rho2;
   double Norm1,Norm2,psf1,psf2;
   double correct2_1,correct2_2; 
   double Etot = fpar[0];
   double X0 = fpar[1];
   double Y0 = fpar[2];
   double sigma1 = fpar[3];
   double sigma2 = fpar[4];
   double ratio = fpar[5];

   double alfa = atan2(Y0,X0);

   eventGetLgdXY(row,col,&x,&y);

   rho2 = SQR(x-X0) + SQR(y-Y0);
   Norm1 = ratio*L2*Etot/(2*3.14159*sigma1*sigma1);
   correct2_1 = L2*(rho2-2.*sigma1*sigma1)/(24.*sigma1*sigma1);
   psf1= Norm1*exp(-0.5*rho2/(sigma1*sigma1))*(1.+correct2_1);

   Norm2 = (1.-ratio)*L2*Etot/(2*3.14159*sigma2*sigma2); 
   correct2_2 = L2*(rho2-2.*sigma2*sigma2)/(24.*sigma2*sigma2);
   psf2= Norm2*exp(-0.5*rho2/(sigma2*sigma2))*(1.+correct2_2);
   return psf1+psf2;
}

