/* public function members (c interface) */

double psfGetE(int row, int col);
double psfGetSigma2(int row, int col); 
int psfSetPars(int npar,double *pars);

double psfGetEI(int row,int col);
double psf1Gss0(int row,int col);
double psf1Gss2(int row, int col);

/*  public function members (fortran interface) */
double psffunc_(int row, int col);
double psfsigma2_(int row, int col);

