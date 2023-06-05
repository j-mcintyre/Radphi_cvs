// CHEES header file
// header file for Keeper Class

#include <stdio.h>

#include "Box.h"
#include <Four.h>
#include <Three.h>
#include <SqMatrix3.h>
#include "fileStructs.h"
extern "C" {
#include <mc_param.h>
}

#ifndef _KEEPER_CLASS_INCLUDED
#define _KEEPER_CLASS_INCLUDED

const int max_kept = 100; /* max # of particles to track */

const int GENERATED_HISTOS = 0;     /* these make accepted plot histo id */
const int ACCEPTED_HISTOS = 1000;   /* equal to generated histo id + 1000 */

class Particle;

class Keeper {

 public:
  mc_codes_t *mCutCodes;   /* link to user input data */
  Particle *mData[max_kept];    /* array of kept particles */
  int mNumParts;                /* number of particles kept */
  
                 Keeper(void);
  void           keep(Particle *);
  unsigned int 	 e852cuts(void);
  
  void           TransformAll(Four *);
  void			 RotateGJtoCM(void);
  void           ConstructRotMatrix(double,double,double,double);
  double	 	 GetTotalEnergy(void);
  
  // input/output and such

  ostream &	print_data(ostream &);
  void 		print_data(char *);
  void		cprint(FILE *, unsigned int); 

  void          InitHistos(void);
  void          FillHistos(int);
  void          WriteHistos(void);
  void          WriteHeader(char *);
  void          SmearE(int);
  void          SmearPosition(position *, int);
  void          mySmear(position*, int);

 private:

  unsigned int  mAcceptBits;

  int 		CircCut(double, double, double, point *);
  void 		getPos(position *, int);

  SqMatrix3 mRotMatrix;
};

#endif


