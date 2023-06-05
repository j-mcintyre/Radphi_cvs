/*-----------------------------------------------------------------------------
Chees.h
 - This is the header file which should be included when using any of the 
   Chees subroutines
 - It may be included in both C and C++ projects
-----------------------------------------------------------------------------*/

/* #ifdef __cplusplus
	extern "C"
	{
          #include "bnlgamsParam.h"
	}
#else
	  #include "bnlgamsParam.h"
#endif */
#ifdef __cplusplus
	extern "C"
	{
          #include "mc_param.h"
	}
#else
	  #include "mc_param.h"
#endif

#include <ntypes.h> 

#ifndef CHEES_INCLUDED
#define CHEES_INCLUDED

#define CHEES_BEAM 0x1
#define CHEES_TARGET 0x2
#define CHEES_RECOIL 0x4
#define CHEES_PRODUCT 0x8
#define CHEES_UNSTABLE 0x10
#define CHEES_FSP 0x20  /* final state particle */
#define CHEES_ESRPART 0x40
#define CHEES_SMEARED 0x100


/*Chees data type definition */


typedef struct
{
  vector4_t p;
  Particle_t id;   /* As defined in particleType.h */
  unsigned int acceptBits;
  int type;  /* see the defines above */
} CheesEntry;



typedef struct
{
  int nVectors;
  int nGammas;
  unsigned int globalAcceptBits;
  vector3_t principleVertex;
  float missingMassSq;
  float genProductMass;
  CheesEntry data[40];
  vector4_t pRawGammas[10];
} cheesInfo_t;



/* function declarations */

  /* main stuff */

void CutChees(char *, mc_param_t *, particles_db_t **);  /* init       */
int  GenerateOne(cheesInfo_t *);                              /* do 1 event */
void ScrubChees(void);                                        /* cleanup    */

  /* utility/support functions */

void ChangeTHisto(int,char*);    /* to change t-distribution histogram */
void setCheesRunNo(int runNo);

#endif
