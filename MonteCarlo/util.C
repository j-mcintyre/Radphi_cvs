/*
   util.C

*/

/*-----------------------------------------------------------------------------
This file contians some utility functions for Chees.

 - createCDF takes a hist id and a filename as input then books and fills a
   histogram with the same limits, bin size, and hist id which is an
   integrated CDF, it returns the header information of the CDF.

 - unfRand returns a random number between 0 and 1. It returns the number in 
   the form of a double.

 - gaussian takes a sigma, a mean, and an x value and returns the value on a 
   gaussian always less than one with sigma and mean which corresponds to that 
   x value

 - GetGaussianDist returns a number after you supply it with a mean
   and a sigma.  This number is picked randomly using a gaussian 
   probability distribution (defined by the sigma and mean you supplied).

 - breitWigner does the same as gaussian but with a Breit-Wigner distribution

 - TDistThree takes a double and a pointer to a mc_transfer_t structure,
   it returns the value of the function defined as tdist #3 (see Rob's input
   parameter specifications) which corresponds to the value of t given as 
   an argument.  Note that this is not neccesarily less than one.

 - doNumIntegration takes a pointer to a function and a pointer to a double.  
   (This function must return a double, takes a double as its first argument, 
   and has a variable argument list.  It cannot belong to a class.)
   It integrates the function until it converges, returns the value of the 
   integration and stores the upperLimit of the integration in the 
   corresponding variable.

 - GetBremDist is just like GetGaussianDist 'cept the probablity distribution
   is 1/x.

 - runEfud is a routine which does a system call to efud.  It takes a the input
   filename to the command as an argument.  If the file is a .bbook file, it 
   truncates the .bbook off the end and puts a .hbook on it, otherwise it just
   addes .hbook to whatever filename you give it.
-----------------------------------------------------------------------------*/

#include "util.h"
#ifdef SOLARIS
#include <ieeefp.h>
extern "C"{
long random();
};
#endif

//////////////////////////////////////////////////////////

wabbitHeader_t createCDF(int histNo, char *filename)
{
  int i;
  float integratedSum = 0.0;
  float runningTotal = 0.0;
  float binWidth = 0.0;
  
  char bbookTitle[MAXSTRING];

  wabbitHeader_t header[1];
  wabbitData_t *data = NULL;

  printf("Integrating histogram %d for CDF...", histNo);

  if(bbGetHisto(histNo, filename, header, &data))
    exit(1);

  sprintf(bbookTitle, "Integrated CDF for histogram %s", header->name);

  bbook1(histNo, bbookTitle, header->nBins1, 
         header->lowerLimit1, header->upperLimit1, header->max);


  /* first get the total sum */
  for(i = 0; i < data->nBins; ++i)
    integratedSum += data->data[i];
  
  /* get more parameters for binning */
  binWidth = (header->upperLimit1 - header->lowerLimit1)/
    header->nBins1;

  /* then fill the CDF */
  for(i = 0; i < data->nBins; ++i)
    {
      runningTotal += data->data[i]/integratedSum;
      bf1(histNo, header->lowerLimit1 + binWidth*i + binWidth/2, runningTotal);
    }

  if(data)
    free(data);

  printf("done\n\n");

  return *header;
}

///////////////////////////////////////////////////////////////

double unfRand(void)
{
  long mrRan;
  double num;
  double f;
  
  mrRan = random();
  num = mrRan; // casting might be dangerous
  f = scalb(num,-31);
  return f;
}

///////////////////////////////////////////////////////////////

double gaussian(double sigma, double mean, double xValue)
{
  return exp(-SQR(xValue-mean)/(2*SQR(sigma)));
}

//////////////////////////////////////////////////////////////

double GetGaussianDist(double mean,double sigma)
{
 /* 
    returns a random number distributed according to a gaussian centered
    on "mean" with width "sigma".  number is always between +-3 sigma.
 */

double x;
int donkeysFly = 1;
while (donkeysFly) {
  x = (mean - 3*sigma) + unfRand()*(6*sigma);

  double testNo = gaussian(sigma,mean,x);

  if ( unfRand() < testNo ) break; 

};

return x;

}

/////////////////////////////////////////////////////////////

double breitWigner(double gamma, double mean, double xValue)
{
  return (SQR(gamma/2)/(SQR(xValue-mean) + SQR(gamma/2)));
}

/////////////////////////////////////////////////////////////

double integrate(double (*func)(double, ...), double *upperLimit)
{
  double xLow = 0.0;
  const double stepSize = 0.01;
  double xHi = stepSize+xLow;
  double area = 0.0;


  double yHi, yLow, lastArea, convergence;

  do
    {
      yHi = (*func)(xHi);
      yLow = (*func)(xLow);
      
      lastArea = area;
      area += stepSize*MIN(yHi, yLow) + 0.5*stepSize*ABS(yHi-yLow);

      convergence = ABS(lastArea-area)/lastArea;
      if(!finite(convergence))
	convergence = 50.0;

      xHi += stepSize;
      xLow += stepSize;
    }while(convergence > 0.000001);

  *upperLimit = xHi-stepSize;

  return area;
}

////////////////////////////////////////////////////////////////////

double GetBremDist(double low,double hi)
{
 /* 
    returns a random number distributed according to 1/x 
    number is always between 'low' and 'hi'.
 */

double x;
int done = 0;
double normalizer = low;
  
  
  do 
  {
    x = low + unfRand()*(hi-low);

    double testNo = normalizer / x ;

    if ( unfRand() < testNo ) done = 1; 

  } while(!done);

  return x;

} /* end getBremDist */


////////////////////////////////////////////////////////////////////

void runEfud(char *outfile)
{
  char systemCall[MAXSTRING];
  char outfileStem[MAXSTRING];

  int isADotBbookFile = FALSE;
  int i,j;

  for(i = 0; outfile[i] != '\0'; ++i)
    if(outfile[i] == '.')
      {
        for(j = 0; j < i; ++j)
          outfileStem[j] = outfile[j];
        outfileStem[i] = '\0';
        
        isADotBbookFile = TRUE;
        break;
      }

  if(isADotBbookFile)
    sprintf(systemCall, "efud -e -i%s -o%s.hbook", outfile, outfileStem);
  else
    sprintf(systemCall, "efud -e -i%s -o%s.hbook", outfile, outfile);

  if(system(systemCall) == -1)
    {
      fprintf(stderr, "Efud has failed to convert %s to hbook format.\n",
              outfile);
      fprintf(stderr, "Try again once program has exited.\n\n");
    }
}
