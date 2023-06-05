#include <stdio.h>
#include <math.h>
#include <rstats.h>
/*
  $Log: average.c,v $
  Revision 1.2  1999/07/05 12:41:42  radphi
  Changes made by gevmage@jlabs1
  updated for 1999

  Revision 1.1  1996/11/21 12:30:32  lfcrob
  Initial revision

*/
static const char rcsid[] = "$Id: average.c,v 1.2 1999/07/05 12:41:42 radphi Exp $";

#define SQR(x) ((x)*(x))

double sum_av[MAXBIN];
double sum_of_squares_av[MAXBIN];
double  npoints_av[MAXBIN];

int zero_av(int bin)

{
  if((bin>=MAXBIN)|(bin<0))
    printf("Invalid bin number %d passed to zero_average\n",bin);
  else{
    sum_av[bin]=0;
    sum_of_squares_av[bin]=0;
    npoints_av[bin]=0;
  }
  return 0;
}

int accumulate_av(int bin,float data)


{
  if((bin>=MAXBIN)|(bin<0))
    printf("Invalid bin number %d passed to accumulate_av\n",bin);
  else{
    sum_av[bin] = sum_av[bin] + data;
    sum_of_squares_av[bin] = sum_of_squares_av[bin] + SQR(data);
    npoints_av[bin] = npoints_av[bin] + 1;
  }
  return 0;
}

int drop_from_av(int bin,float data)


{
  if((bin>=MAXBIN)|(bin<0))
    printf("Invalid bin number %d passed to accumulate_av\n",bin);
  else{
    sum_av[bin] = sum_av[bin] - data;
    sum_of_squares_av[bin] = sum_of_squares_av[bin] - SQR(data);
    npoints_av[bin] = npoints_av[bin] - 1;
  }
  return 0;
}

int get_av(int bin,average_t *average)
{
  float arg;
  if((bin>=MAXBIN)|(bin<0))
    printf("Invalid bin number %d passed to get_average\n",bin);
  else if(npoints_av[bin]>=1){
    average->average = sum_av[bin]/npoints_av[bin];
    arg=sum_of_squares_av[bin]/npoints_av[bin] - SQR(average->average);
    if(arg>0){
      average->sigma = sqrt(arg);
      average->error = average->sigma/sqrt(npoints_av[bin]);
    }
    else{
      average->sigma = 0;
      average->error = 0;
    }
  }
  else
    printf("Not enough data for a fit, bin %d\n",bin);
  return 0;
}
    




