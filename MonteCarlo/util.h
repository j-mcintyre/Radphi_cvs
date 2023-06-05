/* CHEES header file */
/*-----------------------------------------------------------------------------
This file is the header for util.h and any other file that wants to use these 
utilities
-----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include <wabbit.h>
#include <mc_param.h>
}
#include <math.h>
#include <stdarg.h>


#ifndef _UTIL_INCLUDED_
#define _UTIL_INCLUDED_

#define SQR(x) ((x)*(x))                            /* macro for squaring */
#define CUBE(x) ((x)*(x)*(x))                       /* macro for cubing */
#define ABS(x) (((x) > (0.0)) ? (x) : ((-1.0)*(x))) /* absolute value */
#define MIN(x, y) (((x) < (y)) ? (x) : (y))         /* macro for finding min */
#define MAX(x, y) (((x) > (y)) ? (x) : (y))         /* macro for finding max */

/* some constants that everyone can use */
#define PI 3.1415926535897932385
#define E 2.7182818284590452345

#define MAXSTRING 1024
#define TRUE 1
#define FALSE 0

wabbitHeader_t createCDF(int, char *);

double unfRand(void);

double gaussian(double, double, double);
double GetGaussianDist(double ,double );

double GetBremDist(double,double);

double breitWigner(double, double, double);


double integrate(double (*)(double, ...), double *);


void runEfud(char*);

#endif


