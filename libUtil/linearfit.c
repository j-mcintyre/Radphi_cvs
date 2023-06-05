
#include <stdio.h>
#include <math.h>
#include <rstats.h>
/*
  $Log: linearfit.c,v $
  Revision 1.1  1996/11/21 12:30:21  lfcrob
  Initial revision

*/
static const char rcsid[] = "$Id: linearfit.c,v 1.1 1996/11/21 12:30:21 lfcrob Exp $";

double sum_weight[MAXFIT];
double sum_x[MAXFIT];
double sum_y[MAXFIT];
double sum_xsquare[MAXFIT];
double sum_ysquare[MAXFIT];
double sum_xy[MAXFIT];
int npoints_fit[MAXFIT];


int zero_fit(int bin)
{
  if((bin>MAXFIT)|(bin<0)){
    printf("Error in zero_fit, illegal bin %d\n",bin);
    return(0);
  }
  else {
    sum_weight[bin] = 0;
    sum_x[bin] = 0;
    sum_y[bin] = 0;
    sum_xy[bin] = 0;
    sum_xsquare[bin] = 0;
    sum_ysquare[bin] = 0;
    npoints_fit[bin] = 0;
  }
  return(1);
}

int accumulate_fit(int bin,float x,float y,float weight)
{
  if((bin>MAXFIT)|(bin<0)){
    printf("Error in accumulate_fit,Illegal bin %d\n",bin);
    return(0);
  }
  else if(weight<=0){
    printf("Error in accumulate_fit,Illegal weight %f\n",bin);
    return(0);
  }
  else {
    npoints_fit[bin]++;
    sum_weight[bin] = sum_weight[bin] + (double) weight;
    sum_x[bin] = sum_x[bin] + ((double)x)*((double)weight);
    sum_y[bin] = sum_y[bin] + ((double)y)*((double)weight);
    sum_xy[bin] = sum_xy[bin] + ((double)x)*((double)y)*((double)weight);
    sum_xsquare[bin] = sum_xsquare[bin] + pow((double)x,2.0)*((double)weight);
    sum_ysquare[bin] = sum_ysquare[bin] + pow((double)y,2.0)*((double)weight);
  }
  return(1);
}
    
int get_fit(int bin,linearfit_t *fit)

{
  double det,dslope,dintercept;

  if((bin>MAXFIT)|(bin<0)){
    printf("Error in get_fit,Illegal bin %d\n",bin);
    return(0);
  }
  else if(npoints_fit[bin]<2){
    printf("Error in get_fit,Not enough points for a fit in bin %d\n",bin);
    return(0);
  }
  else{
    det = sum_xsquare[bin]*sum_weight[bin]-pow(sum_x[bin],2.0);
    dslope = (sum_weight[bin]*sum_xy[bin] - sum_x[bin]*sum_y[bin])/det;
    fit->slope = (float) dslope;
    dintercept = (sum_xsquare[bin]*sum_y[bin] - sum_x[bin]*sum_xy[bin])/det;
    fit->intercept = (float) dintercept;
    if(sum_weight[bin]>2)
      fit->chisquare = (sum_ysquare[bin] + (float) pow((double)fit->slope,2.0)*sum_xsquare[bin] - 2*(fit->slope)*
	sum_xy[bin] -2*(fit->intercept)*sum_y[bin] + 2*(fit->slope)*(fit->intercept)*sum_x[bin]
	  + (float) pow((double)fit->intercept,2.0)*sum_weight[bin])/(sum_weight[bin]-2.0);

/*	  {
	  	fit->chisquare = (float)
          (sum_ysquare[bin] 
      	  + pow(dslope,2.0)*sum_xsquare[bin] 
      	  - 2*(dslope)*sum_xy[bin] 
      	  - 2*(dintercept)*sum_y[bin] 
      	  + 2*(dslope)*(dintercept)*sum_x[bin]
	  	  + pow(dintercept,2.0)*sum_weight[bin])
         / (sum_weight[bin]-2.0);
         printf("\nchi2(%d)=%f",bin,fit->chisquare);
	 }*/
    else
      {
      	fit->chisquare=0;
      }
    fit->slope_err = sqrt(sum_weight[bin]/det);
    fit->intercept_err = sqrt(sum_xsquare[bin]/det);
  }
  return(1);
}

