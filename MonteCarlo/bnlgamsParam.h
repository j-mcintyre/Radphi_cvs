/*
   bnlgamsParam.h 5.2 
   7/27/97 18:52:19
   /home/lemond/e852/source/include/SCCS/s.bnlgamsParam.h
   @(#)bnlgamsParam.h	5.2
*/

#include <particleType.h>

#ifndef BNLGAMSPARAM_H_INCLUDED
#define BNLGAMSPARAM_H_INCLUDED


#define BEAM_MOMENTUM_FIXED 0x1
#define BEAM_MOMENTUM_GAUSSIAN 0x2
#define BEAM_MOMENTUM_FLAT 0x4
#define BEAM_MOMENTUM_BREM 0x8
#define BEAM_MOMENTUM_BNL_94 0x10

#define BEAM_SPATIAL_NONE 0x1
#define BEAM_SPATIAL_GAUSSIAN 0x2
#define BEAM_SPATIAL_CIRCLE 0x4

#define BEAM_ANGULAR_NONE 0x1
#define BEAM_ANGULAR_DIVERG 0x2

#define MASS_FLAT 0x1
#define MASS_BW 0x2
#define MASS_GAUSSIAN 0x4
#define MASS_DELTA 0x8
#define MASS_USER 0x10

#define DECAY_THETA_YLM 0x1
#define DECAY_THETA_FLAT 0x2
#define DECAY_THETA_UNIFORM 0x4
#define DECAY_THETA_DALITZ 0x8
#define DECAY_THETA_USER 0x10
#define DECAY_THETA_CORRELATED 0x20

#define DECAY_THETA_MASK 0xff

#define DECAY_PHI_FLAT 0x100
#define DECAY_PHI_UNIFORM 0x200
#define DECAY_PHI_USER 0x400

#define DECAY_PHI_MASK 0xff00


#define TDIST_1 0x1
#define TDIST_2 0x2
#define TDIST_3 0x4
#define TDIST_USER 0x8
#define TDIST_DELTA 0x10
#define TDIST_FLAT 0x20

#define GEOM_E852 0x1
#define GEOM_CEBAF 0x2
#define GEOM_GAMS2000 0x4
#define GEOM_GAMS4000 0x8

#define OUT_NONE 0x1
#define OUT_C 0x2
#define OUT_CPLUSPLUS 0x4
#define OUT_FORTRAN 0x8
#define OUT_ASCII 0x10


typedef struct{
  char filename[1024];
  int id;
} histogramDef_t; 

typedef struct{
  int pCode;
  float p;
  float bremThetaMax;
  float pWidth;
  int spreadCode;
  float sigma[2];  
  float radius;  
  int angularCode;
  float diverg[2];
} bnlgams_beam_t;


typedef struct{
  int transferCode;  
  float exponent;
  float exchangeMass;
  float spinPower;
  float tPower;
  float t;
  float tWidth;
  histogramDef_t histo;
} bnlgams_transfer_t; 

typedef struct{
  int nProducts;
  float branchingRatio;
  int product[4];
} decay_t;



typedef struct{
  int code;
  histogramDef_t thetaHisto;
  histogramDef_t phiHisto;
  float theta;
  float thetaWidth;
  float phi;
  float phiWidth;
} spinDecay_t;

typedef struct{
  char name[16]; 
  Particle_t id;
  int isPhoton;
  int charge;
  float lifetime;
  int massCode;
  float mass;
  float sigma;
  histogramDef_t massHisto;
  int spin;
  spinDecay_t twoBody;
  spinDecay_t threeBody;
  int nDecays;
  int returnData;
  decay_t decay[10];
} particle_db_t;

typedef struct{
  int nParticles;
  particle_db_t particle[1];
} particles_db_t;

typedef struct{
  int nIterations;
  int doHisto;
  int geometry;
  int returnAll;
  int outFileType;
  char outFile[1024];
  int doEnergySmear;
  int doPositionSmear;
  float separationCut;
  float energyCut;
  float energySmearFloor;
  float energySmearStat;
  float positionSmear[2];
  int fermiTarget;
} bnlgams_codes_t;

typedef struct{
  int beam;
  int target;
  int recoil;
  int product;
} bnlgams_reaction_t;

typedef struct{
  bnlgams_reaction_t reaction;
  bnlgams_beam_t beam;
  bnlgams_transfer_t transfer;
  bnlgams_codes_t codes;
} bnlgams_param_t;

#define PARAM_OK 0
#define PARAM_ERROR 1
#define LOAD_PARTICLE_FAIL -1
#define BOOLEAN_FAIL -1

#define PARTICLE_DB "particles.xrdb"

int  cheesParamAddFile(char *);
int initParam(char *);
int loadParam(bnlgams_param_t *param,particles_db_t **db);
void flushParam(void);

#endif
