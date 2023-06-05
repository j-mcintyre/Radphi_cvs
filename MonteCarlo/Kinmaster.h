// CHEES header file
/* Kinmaster object header

Kinmaster is in charge of kinematic quantities that are not associated
with the Particle class.

PCF 22 June 1995

*/

extern "C" {
#include <mc_param.h>
}

#include "Particle.h"
#include "NosPicker.h"

#ifndef _KINMASTER_CLASS_INCLUDED_
#define _KINMASTER_CLASS_INCLUDED_

#define TDIST3_HIST 30
#define BREMTHETA_HIST 31
class Kinmaster
{
  public:
  
  // data //
    
    /* 
      These are links to the beam and momentum transfer (t) info,
      which was read in from the input file into global memory.
    */
    
    
  mc_beam_t     *mBeam;
  mc_transfer_t *mTransfer;
  
  double mTotCMEnergy;  /* total energy available in the overall CM */
  double mCurrentT;     /* current momentum transfer for this reaction */
  double mThetaC;       /* used in brem beam angular distribution */
  
  NosPicker *mUserTPicker;      /* for user specified t */
  NosPicker *mUserFermiPicker;  /* for user specified Fermi momentum
                                   distribution */
  NosPicker *mUserBremThetaPicker; /* for brem beam divergence */
  // member functions //

  Kinmaster();              // constructor
  ~Kinmaster();             // destructor

  // InitKinmaster does some initialization bullshit needed for choosing t, 
  // and it is stuff that can't be put into the constructor because it needs
  // mTransfer to be set before it can do it.  It also initializes the 
  // NosPicker
  
  void InitKinmaster();

  
  void 	   FillCMVectors(double, double, double, Particle*, Particle*);
  void     Fill2BodyGJInverse(double, double, double, Particle*, Particle*);
  double   ChooseT(double,double,Particle*,Particle*);
  double   GetProductCMEnergy(double,Particle*,Particle*);

  void     PrepareBeam(Particle*); //  it fills beam momenta
  void     PrepareTarget(Particle*); //  it fills target momenta

  void LorentzXForm1(Four *, Four *);    // lorentz transforms one particle

private:

    // data //
  
  double     mTMin;       /* minimum allowable t, computed in InitChooseT */
  double     mTMax;       /* maximum allowable t, computed in InitChooseT */

    // member functions //

  // The following functions take an xValue as an argument and return the 
  // corresponding value on a t-distribution.  TDistOne and TDistTwo are 
  // constructed so that their maximum value is 1.  TDistThree is not done
  // this way (too complicated).
  
  inline double TDistOne(double);       // exp{-b*t}
  inline double TDistTwo(double);       // t*exp{-b*t}
  inline double TDistThree(double);     // t^n*(t-tMin)^w*e^(-b*t)/(t+mX^2)^2
  inline double BremThetaFunc(double);   // x / (thetaC^2 + x^2)^2

  double DoNumIntegration(void);
 
  void FillTDistThree(double, double);  // fills a histogram
  void FillBremDist(double,double);  // fills brem beam theta histo
  inline double ChooseGuessT(void);     // chooses a first guess within a range

  inline double GetZVertex(void);            // choose vertex Z coordinate 
  inline void GetTargetCenter(dvector3_t*);   // Get targ cent. from targets
                                             // header 
                                             
  
};

#endif




