
#ifndef RSTATS_H_INCLUDED
#define RSTATS_H_INCLUDED

static const char rcsid_rstats_h[]="$Id: rstats.h,v 1.1 1996/11/21 12:31:14 lfcrob Exp $";


#define MAXFIT 100000
#define MAXBIN 100000

typedef struct{
  float slope;
  float intercept;
  float slope_err;
  float intercept_err;
  float chisquare;
}linearfit_t;

typedef struct{
  float average;
  float error;
  float sigma;
}average_t;

int accumulate_av(int,float);
int zero_av(int);
int drop_from_av(int,float);
int get_av(int,average_t *);

int zero_fit(int);
int accumulate_fit(int bin,float x,float y,float weight);
int get_fit(int bin,linearfit_t *fit);

#endif
