/*
 * clib.h
 *
*/


#ifndef clibH
#define clibH

#include <stdio.h>
#include <time.h>
#include <ntypes.h>
#include <stdarg.h>

#define TIME_SCALE 250000
#define TWO_DAYS 172800
#define MAXLINE 136
#define MAXCNT   30
#define ENDLINE  77 

/* Misc stuff */

char *Ctime(const time_t  *);


/* Byte order manipulation */

void flip(unsigned char *a, unsigned char *b);
void flip_integer(unsigned char *integer);
void flip_short(unsigned char *short_int);
void flipint(char *c);
void flipshort(char *c);

/* Gaussian distributions */

void GaussSet(float mean,float width);
float gasdev(void);
float gauss(float mean,float width);
void gaussRndm2(double mean1,double mean2,double sigma1,double sigma2,double *rnd1,double *rnd2);

/* Misc numeric stuff */

int is_NaN(double);


/* Vector manipulation */

float v3mag(const vector3_t* v3);
float v3magsq(const vector3_t* v3);

float v4mag(const vector4_t* v4);
float v4magsq(const vector4_t* v4);

vector3_t v3sub(const vector3_t *v1,const vector3_t *v2);
vector4_t v4sub(const vector4_t *v1,const vector4_t *v2);

void v3sum(const vector3_t *v1,const vector3_t *v2,vector3_t *sum);
void v4sum(const vector4_t *v1,const vector4_t *v2,vector4_t *sum);

void add3(vector3_t *v,const vector3_t *plus);
void add4(vector4_t *v,const vector4_t *plus);

vector3_t v3add(const vector3_t *v1,const vector3_t *v2);
vector4_t v4add(const vector4_t *v1,const vector4_t *v2);

void v3diff(const vector3_t *v1,const vector3_t *v2,vector3_t *diff);
void v4diff(const vector4_t *v1,const vector4_t *v2,vector4_t *diff);

void sum4vectors(int n,const vector4_t *fourvecs[],vector4_t *sum);
void sum3vectors(int n, vector3_t *threevecs[],vector3_t *sum);

/*
  Use n4VectSum as vector4_t pTot = n4VectSum(int n,vector4_t *p1,
                                    vector4_t *p2,....,vector4_t *pn);
*/
vector4_t BetaGamma(const vector4_t *p);

vector4_t n4VectSum(int nVect,...);


vector3_t negative(const vector3_t *v);
vector4_t scalerMult4V(vector4_t *vect, float scale);

void p3vector(const vector3_t *v,FILE* fp);
void p4vector(const vector4_t *v,FILE* fp);

void p4vectors(FILE* fp,int npart,const vector4_t *vectors[],const char* title);

double DotProduct3(const vector3_t *p1,const vector3_t *p2);
vector3_t CrossProduct3(const vector3_t *p1,const vector3_t *p2);

void make4vector(const vector3_t *v3,float mass,vector4_t *v4);

/*
  Use nPartMass as float mass = nPartMass(int n,vector4_t *p1,
                                vector4_t *p2,....,vector4_t *pn);
*/

float nPartMass(int nPart,...); /* the ... is a list of arguments, each one
				   a pointer to a vector4_t */

float EffMass(int npart,const vector4_t *fourvecs[]);

/*
  effMass() Returns "signed" mass (i.e. for MMSQ)
  */
float effMass(int npart,const vector4_t *fourvecs[]);
vector4_t lorentz(const vector4_t *beta,const vector4_t *pin);

float ctheta(const vector3_t*vector1,const vector3_t*vector2); /* cos(angle) between v1 and v2 */

void ScalerMult(vector3_t *vout,float s,const vector3_t *vin);

#endif
/* end file */
