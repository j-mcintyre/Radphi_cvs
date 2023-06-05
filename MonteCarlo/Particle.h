// CHEES header file
// header for Particle class
// PCF, 16 June 1995
/////////////////////////////////////////////////////

#include <iostream>
using namespace std;

extern "C" {
#include <mc_param.h>
}
#include "Four.h"
#include "NosPicker.h"
#include "ReadUserIntensity.h"

#ifndef _PARTICLE_CLASS_INCLUDED_
#define _PARTICLE_CLASS_INCLUDED_

#define HIST_Y10 10
#define HIST_Y11 11
#define HIST_Y20 20
#define HIST_Y21 21
#define HIST_Y22 22
#define HIST_Y40 40


const int max_daughters = 3;  /* max # of daughter particles */

//////////////////////////////////////////////////////

class Particle
{

public:

/*--------------------------- data ----------- */
  particle_db_t *mInfo;  /* 
                              Link to particle info; lifetime,
                              name, charge, decay info, you name it. 
                              This info is read from the input file
                              so we just give the Particle a pointer to it.
                         */  

  dvector3_t mOrigin;           /* where the Particle was born */
  
  double mMass;     	   /* the particle's mass */ 
  Four mFour;     	   /* the particle's four vector */
  
  Particle *mDaughter[max_daughters];  /* pointers to daughter particles */
  int mNumDaughters;             /* how many daughter particles */
  
  Particle *mParent;

  /*  double cosThetaLab;     cos theta of momentum vector in lab frame
  double phiLab;              phi of momentum vector in lab frame  */
  
  double mCosTheta;        /* cosTheta of decay in G-J frame */
  double mPhi;             /* phi of decay in G-J frame */
  
  unsigned int	mPersonalAcceptBits;  /* flags for cuts 0 = good, 1 = bad */
                                      /* see "cutmasks.h" */

/*--------------------------- member functions */
  
  Particle();                	  // constructor
  ~Particle();                    // destructor
  void InitParticle(void);        // basically just inits the NosPickers

  double Mass(void);              // easy way to calculate sqrt(E^2 - p^2)
  double MassSq(void);            // get the mass of the particle squared
  double Momentum(void);          // get the magnitude of the three-momentum
  void Decay(Keeper *, 
             Keeper *,
             int      );          // fill in the Particle daughter four vectors
  void SetEnergy(double); 	  // sets the particle's energy
  double Energy(void); 		  // gets particle's energy
  
  void ChooseMass(void);          // selects a mass according to specifications
  void operator=(Particle &); 	  // assigns values in one Particle to another
  friend ostream &operator<<(ostream &, Particle &);   // prints out a Particle
  
private:

/*----------------- private data --------------*/


   NosPicker *mUserMassPicker;
   NosPicker *mUserCosThetaPicker;
   NosPicker *mUserPhiPicker;
  
   // stuff for intensity plots using correlated m,costheta,phi distributions 
   float *mI;
   int mNMassBins;
   int mNCthBins;
   int mNPhiBins;
   float mMassLow;
   float mMassBinSize;
   float *mScale;
   //

   void decay2(int); 	              // decays 2 daughter particles
   void decay2_select(int);         
     /* decay2_select is Craig Steffen's decay function that decays the
	particle, but does statistical selections of other things. 
	Implemented 22JAN99     - Craig Steffen  */
   void decay3(void); 		      // decays 3 dauther particles
   inline double ChooseGuessMass();   // chooses a first guess within a range
   
};

#endif














