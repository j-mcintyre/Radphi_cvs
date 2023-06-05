/* 

Kinmaster object 

Kinmaster is in charge of kinematic quantities that are not associated
with the Particle class.

PCF 22 June 1995

*/

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#ifdef SOLARIS
#include <ieeefp.h>
#endif

#include "Kinmaster.h"
#include <mc_param.h>
#include "util.h"
#include "targets.h"
#include <SqMatrix3.h>
#include <Three.h>

extern mc_param_t *gInParam;
extern Particle *gBeam;
extern Particle *gTarget;
extern Particle *gRecoil;

char *original = "newPain.bbook";
char *cdf = "joy.bbook";

char *brem_dist_file = "bremdist.bbook";
char *brem_CDF_file  = "bremcdf.bbook";

#define FERMI_FILE "/fermiCDF.bbook"
static char *fermiFile;

//////////////////////////////////////////////////////////////////

Kinmaster::Kinmaster()  //constructor
{
   mTotCMEnergy = 0.0;

  /* set up the NosPickers, which choose random values
     according to user specified probability distributions */

   /* allocate space for the t distribution NosPicker */
   mUserTPicker = NULL;
   mUserTPicker = new(NosPicker);
   assert(mUserTPicker != NULL);
   
   /* allocate space for fermi momentum NosPicker */
   mUserFermiPicker = NULL;
   mUserFermiPicker = new(NosPicker);
   assert(mUserTPicker != NULL);
   
    /* allocate space for brem theta  NosPicker */
   mUserBremThetaPicker = NULL;
   mUserBremThetaPicker = new(NosPicker);
   assert(mUserBremThetaPicker != NULL);


}


//////////////////////////////////////////////////////////////////

Kinmaster::~Kinmaster()  //destructor
{
  // clean up droppings if user has specified TDIST_3 or brem beam

  char systemCall[MAXSTRING];

  if(mTransfer->transferCode == TDIST_3)
    {
      sprintf(systemCall, "rm %s %s", original, cdf);
      if(system(systemCall) == -1)
	{
	  fprintf(stderr, "Failed to remove files %s and %s created ");
	  fprintf(stderr, "for construction of t-dist #3\n");
	  fprintf(stderr, "Try again once program has exited\n\n");
        }
    }

  if(mBeam->pCode == BEAM_MOMENTUM_BREM)
    {
      sprintf(systemCall, "rm %s %s", brem_dist_file, brem_CDF_file);
      if(system(systemCall) == -1)
	{
	  fprintf(stderr, "Failed to remove files %s and %s created ");
	  fprintf(stderr, "for construction of brem beam\n");
	  fprintf(stderr, "Try again once program has exited\n\n");
        }
    }

}

//////////////////////////////////////////////////////////////////

void Kinmaster::InitKinmaster()
{
  // One time initialization that has to wait until after
  // the input file is in memory.

  // Fermi momentum: ///////////////////////////////////////
  //
  //
  	if(gInParam->codes.fermiTarget)
  	{
  	  // User wants fermi momentum.  Initialize: 
  	  
         if (!fermiFile) // add histo directory name to bbook filename
         {
          char *dir;
	  char *subdir="/MonteCarlo/";
          dir = getenv("RADPHI_CONFIG");
          assert(dir);
          fermiFile = (char *)calloc(strlen(dir)+strlen(FERMI_FILE)+strlen(subdir)+1,
                        sizeof(char));
          strcpy(fermiFile,dir);
          strcat(fermiFile,subdir);
          strcat(fermiFile,FERMI_FILE);
         };
           
  	  if(mUserFermiPicker->InitNosPicker(1 /* histo # */, fermiFile))
  	  {
  	  	
	  /* error message handling has already been done */
	  fprintf(stderr, "Error reading fermiCDF.bbook (Kinmaster.C)\n");
	  fprintf(stderr, "Exiting...\n");
	  exit(EXIT_FAILURE);
	  };
	  
	  // Now mUserFermiPicker is all set to go.
	
	}; /* end if fermiTarget*/ 	
  //
  //
  /////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////
  // t distribution:
  // If the transferCode calls for a user-specified distribution, then  
  // initialize the histogram.  If the transferCode wants Tdist #3, then it 
  // is too difficult to pick values off it the same way as for Tdist #1 and
  // #2 so we have to integrate the babe numerically, save it to a file, then
  // call the same routines that we use for user-specified functions

  double upperLimit = 0.0;
  double binWidth = 0.01;
  double histoLimit = 0.0;

  /* first set mTMin as if it were zero */
  mTMin = 0.0;
	
  switch(mTransfer->transferCode)
    {
    case TDIST_USER:
      if(mUserTPicker->InitNosPicker(mTransfer->histo.id, mTransfer->histo.filename))
	{
	  /* error message handling has already been done */
	  fprintf(stderr, "Check your input parameter file for for possible ");
	  fprintf(stderr, "sources of error\n");
	  fprintf(stderr, "Exiting...\n\n");
	  exit(EXIT_FAILURE);
	}
      break;
      
    case TDIST_3:
      upperLimit = Kinmaster::DoNumIntegration();
      
      /* now we need to round upperLimit to the nearest bin so that we don't
         get a spike in the middle of the tdistribution */
      if(upperLimit/binWidth - (int)(upperLimit/binWidth) >= 0.5)
	histoLimit = (double)((int)(upperLimit/binWidth)+1)*binWidth;
      else
	histoLimit = (double)((int)(upperLimit/binWidth))*binWidth;
      
      bbook1(TDIST3_HIST, "t-Distribution #3", (int)(histoLimit/binWidth), 
	     0, upperLimit, 0);
      FillTDistThree(histoLimit, binWidth);
      bbsave(TDIST3_HIST, original);
      createCDF(TDIST3_HIST, original);
      bbsave(TDIST3_HIST, cdf);
      if(mUserTPicker->InitNosPicker(TDIST3_HIST, cdf))
	{
	  /* error message handling has already been done */
	  fprintf(stderr, "Creation of CDF for tdist #3 joy.bbook has");
	  fprintf(stderr, "failed\n");
	  fprintf(stderr, "This error should not occur\n");
	  fprintf(stderr, "Contact Chees expert to fix bug\n");
	  fprintf(stderr, "Exiting...\n");
	  exit(EXIT_FAILURE);
	}
      break;
      
    default:
      break;
      
    } /* end of transfer code switch */     
	
	//
	//
	/////////////////////////////////////////////

    /////////////////////////////////////////////
    //
    //  brem beam
    // 
    // the theta distribution for the brem beam is
    // complicated to it pays to make a CDF histo
    // here too (as in tdist 3 above).
	
	if (mBeam->pCode == BEAM_MOMENTUM_BREM)
	{    
      histoLimit = 0.0;
  
	  mThetaC = (0.000511) / (mBeam->p);  // mass of electron over Emax
	 									  // is characteristic brem photon angle
	  
	  										 
      upperLimit = mBeam->bremThetaMax; // it's in units of mThetaC
      							
      							
      upperLimit *= mThetaC; // convert to radians
      
      binWidth = upperLimit/100; // 100 bins
      
      /* now we need to round upperLimit to the nearest bin so that we don't
         get a spike in the middle of the tdistribution */
      
      if(upperLimit/binWidth - (int)(upperLimit/binWidth) >= 0.5)
	    histoLimit = (double)((int)(upperLimit/binWidth)+1)*binWidth;
      else
	    histoLimit = (double)((int)(upperLimit/binWidth))*binWidth;
      
      bbook1(BREMTHETA_HIST, "brem beam theta dist", (int)(histoLimit/binWidth), 
	     0, upperLimit, 0);
      FillBremDist(histoLimit, binWidth);
      bbsave(BREMTHETA_HIST, brem_dist_file);
      createCDF(BREMTHETA_HIST, brem_dist_file);
      bbsave(BREMTHETA_HIST, brem_CDF_file);
      
      if(mUserBremThetaPicker->InitNosPicker(BREMTHETA_HIST, brem_CDF_file))
	  {
	  /* error message handling has already been done */
	  fprintf(stderr, "Creation of CDF for brem beam theta has");
	  fprintf(stderr, "failed\n");
	  fprintf(stderr, "This error should not occur\n");
	  fprintf(stderr, "Contact Chees expert to fix bug\n");
	  fprintf(stderr, "Exiting...\n");
	  exit(EXIT_FAILURE);
	   }
    
    }
    //
    //
    /////// end brem beam init
}
//////////////////////////////////////////////////////////////////

void Kinmaster::PrepareTarget(Particle *theTarget)
{

	    double targetMass = theTarget->Mass(); // we'll use this below.

if(gInParam->codes.fermiTarget)
{
	// If fermi momentum is on, assign the target a momentum according
	// to the distribution in fermiCDF.bbook.
		
	// first pick a value for the magnitude of Pf:
	
		double momentum = mUserFermiPicker->Pick();
	
	// (Hey, that was easy, don't you love C++?)
	// randomly pick an orientation for the momentum
	  
	  double cosTheta,sinTheta,phi;
      do
      {
	    cosTheta = (2.0 * unfRand()) - 1.0;
	    sinTheta = sqrt(1.0 - (cosTheta*cosTheta)); 
        phi = 2.0*PI*unfRand(); 
	  }
       while (!finite(sinTheta)); 

	// now we have magnitude and direction for Pf.  Assign components.
	
	    float px = momentum * sinTheta * cos(phi);
	    float py = momentum * sinTheta * sin(phi);
	    float pz = momentum * cosTheta;
	    	    
	    theTarget->mFour.Set(px,py,pz,targetMass);
	    
	// done.
}
else
{
	// Target always has zero momentum. 
	
	theTarget->mFour.Set(0.0,0.0,0.0,targetMass);
	
};

}

//////////////////////////////////////////////////////////////////

void Kinmaster::PrepareBeam(Particle *theBeam)
{
static double pbeam_x,pbeam_y,pbeam_z,ebeam;
dvector3_t targetCenter;
double theta,origin_x,origin_y,origin_z;


// Set up the beam according to user input specifications.
// Assumes the beam mass is set in mMass.

// If they've chosen BEAM_MOMENTUM_BNL, another routine
// is called to prepare the beam momentum AND position-- so we'll bypass 
// the part of this routine that generates the beam position and momentum:

if ((mBeam->pCode) != BEAM_MOMENTUM_BNL_94)
{
  
  // First find the magnitude of the beam momentum.
  
  double magBeamP;

   switch (mBeam->pCode)
   {
     case BEAM_MOMENTUM_FIXED:
        magBeamP = mBeam->p;
       break;
 
     case BEAM_MOMENTUM_GAUSSIAN:
    	magBeamP = GetGaussianDist(mBeam->p,mBeam->pWidth); 
       break;

     case BEAM_MOMENTUM_FLAT:  /* p - sig < p_z < p + sig */

        magBeamP = (mBeam->p - mBeam->pWidth) + unfRand()*2*(mBeam->pWidth);
       break;

     case BEAM_MOMENTUM_BREM:
        magBeamP = GetBremDist(  0.75 * (mBeam->p) /* min */,
                               0.95 * (mBeam->p) /* max */);
       break;

     default: /* they messed up ! */
       fprintf(stderr, "Invalid beam momentum specification. (Kinmaster::PrepareBeam)\n");
       fprintf(stderr, "Check your input parameter file.\n\n");
       exit(EXIT_FAILURE);
       break;

   } /* end beam momentum switch */


// Now pick x & y momenta for the beam, if the user
// has specified angular beam divergence.
//
// Brem beam implies divergence as well.
//

double theta_x, theta_y;

   switch (mBeam->angularCode)
   {

     case BEAM_ANGULAR_NONE:
        pbeam_x = 0.0;
        pbeam_y = 0.0;
        pbeam_z = magBeamP;
       break;

     case BEAM_ANGULAR_DIVERG:
        theta_x = GetGaussianDist(0.0,mBeam->diverg[0]); 
        theta_y = GetGaussianDist(0.0,mBeam->diverg[1]); 

        // pbeam_z has already been set //

        pbeam_x = magBeamP * sin(theta_x);
        pbeam_y = magBeamP * sin(theta_y);
	    pbeam_z = sqrt( SQR(magBeamP) - SQR(pbeam_x) - SQR(pbeam_y) ); 
       break;

     default:
       if((mBeam->pCode) != BEAM_MOMENTUM_BREM){
	 fprintf(stderr, "Invalid beam angular divergence specification.\n");
	 fprintf(stderr, "Check your input parameter file.\n\n");
	 exit(EXIT_FAILURE);
       }
       break;
   } /* end angular switch */
   
   
 // If they are using a brem beam, handle photon direction independently 
 // independently from the above.
 // Theta,phi refer to the angles w.r.t. z axis
 
    if ((mBeam->pCode) == BEAM_MOMENTUM_BREM)
    {
    
    double bremTheta = 100.0;
    double bremPhi;              // polar angles of photon beam w.r.t. z-axis 
    
    
     // choose theta
     
    double thetaC = (0.000511/mBeam->p);  // electron mass by electron energy
    double thetaMax = thetaC * mBeam->bremThetaMax;  // thetaMax in units of thetaC 
    while (bremTheta > thetaMax)
     { bremTheta = mUserBremThetaPicker->Pick(); }
    
    // choose phi
    
    bremPhi = 2.0*PI*unfRand();  // uniform between 0,2pi
    
    // divide into components
    
    pbeam_z = magBeamP * cos(bremTheta);
    pbeam_x = magBeamP * sin(bremTheta) * cos(bremPhi);
    pbeam_y = magBeamP * sin(bremTheta) * sin(bremPhi);
    
    }
 
   
// Now find where the beam is in space (user specified).
//  These are the vertex x,y coordinates.
// Also choose the vertex Z coordinate (currently not user specified).
// Put this info in mOrigin of the beam.


   this->GetTargetCenter(&targetCenter);
   origin_z = GetZVertex();

   switch (mBeam->spreadCode)
   {

     case BEAM_SPATIAL_NONE:
        origin_x = targetCenter.x;
        origin_y = targetCenter.y; 
       break;

     case BEAM_SPATIAL_GAUSSIAN:
        origin_x = GetGaussianDist(targetCenter.x,mBeam->sigma[0]);
        origin_y = GetGaussianDist(targetCenter.y,mBeam->sigma[1]);
       break;
     
     case BEAM_SPATIAL_CIRCLE:      // randomly populated circle
        theta = 2*PI *  unfRand();  // pick random angle
        origin_x = targetCenter.x + unfRand() * mBeam->radius * cos(theta); 
        origin_y = targetCenter.y + unfRand() * mBeam->radius * sin(theta); 
       break;

     default:
       fprintf(stderr, "Invalid beam spatial spread specification.\n");
       fprintf(stderr, "Check your input parameter file.\n\n");
       exit(EXIT_FAILURE);
       break;

   } /* end spatial switch */



}/* end of if ! BEAM_MOMENTUM_BNL */
else
#ifdef DO_E852
{
    // They've chosen BEAM_MOMENTUM_BNL.
    // Call the routine that prepares the beam.
    

	dvector3_t beamPosition, beamMomentum;
	
	while(E852beamPicker(E852_18_1994,&beamMomentum,&beamPosition));
	
	pbeam_x = beamMomentum.x;
	pbeam_y = beamMomentum.y;
	pbeam_z = beamMomentum.z;
	
	origin_x = beamPosition.x;
	origin_y = beamPosition.y;
	origin_z = beamPosition.z;    

}; /* end of BEAM_MOMENTUM_BNL */
#else
{
fprintf(stderr,"E852 Beam is not supported (Kinemaster.C)\n");
};
#endif

// Now that we've calculated everything, we can actually set
// the four vector

   // e^2 = p^2 + m^2

ebeam = sqrt( 
              SQR(pbeam_x) +
              SQR(pbeam_y) +
              SQR(pbeam_z) +
              SQR(theBeam->mMass)
            );

theBeam->mFour.Set(
                    pbeam_x    /* px */
                   ,pbeam_y    /* py */
                   ,pbeam_z    /* pz */
                   ,ebeam      /* E  */
                  );

theBeam->mOrigin.x = origin_x;  // set the vertex 
theBeam->mOrigin.y = origin_y;
theBeam->mOrigin.z = origin_z;

}

//////////////////////////////////////////////////////////////////

inline void Kinmaster::GetTargetCenter(dvector3_t *center)
{

  switch(gInParam->codes.geometry)
  {
    case GEOM_E852:
#ifdef DO_E852
       center->x = e852_target_center_x;
       center->y = e852_target_center_y;
       center->z = e852_target_center_z;
#else
       fprintf(stderr, "E852 geometry not yet implemented (Kinmaster::GetTargetCenter.\n");
       exit(EXIT_FAILURE);
#endif
      break;

    case GEOM_CEBAF_328:
    case GEOM_CEBAF_620:
       center->x = 0.0;
       center->y = 0.0;
       center->z = 0.0;
      break;
     
    case GEOM_GAMS2000:
       fprintf(stderr, "GAMS geometry not yet implemented (Kinmaster::GetTargetCenter.\n");
       exit(EXIT_FAILURE);
      break;

    case GEOM_GAMS4000:
       fprintf(stderr, "GAMS geometry not yet implemented (Kinmaster::GetTargetCenter).\n");
       exit(EXIT_FAILURE);
      break;

    default: 
       fprintf(stderr, "Invalid geometry specification.\n");
       fprintf(stderr, "Check your input parameter file.\n\n");
       exit(EXIT_FAILURE);
      break;

  }; /* end switch */
}


//////////////////////////////////////////////////////////////////

inline double Kinmaster::GetZVertex(void)
{
double answer;
  switch(gInParam->codes.geometry)
  {
    case GEOM_E852:
#ifdef DO_E852
       answer =  (e852_target_low_z + unfRand()*(e852_target_hi_z -
                                              e852_target_low_z ));
#else
       fprintf(stderr, "E852 geometry not yet implemented (Kinmaster::GetZVertex).\n");
       exit(EXIT_FAILURE);
#endif
      break;

    case GEOM_CEBAF_328:
    case GEOM_CEBAF_620:
       answer = 0.0;
      break;
     
    case GEOM_GAMS2000:
       fprintf(stderr, "GAMS geometry not yet implemented (Kinmaster::GetZVertex).\n");
       exit(EXIT_FAILURE);
      break;

    case GEOM_GAMS4000:
       fprintf(stderr, "GAMS geometry not yet implemented (Kinmaster::GetZVertex).\n");
       exit(EXIT_FAILURE);
      break;

    default: 
       fprintf(stderr, "Invalid geometry specification.\n");
       fprintf(stderr, "Check your input parameter file.\n\n");
       exit(EXIT_FAILURE);
      break;

  }; /* end switch */

 return answer;
}

//////////////////////////////////////////////////////////////////

void Kinmaster::Fill2BodyGJInverse(  
                           double inTotalEnergy
                          ,double inCosTheta
                          ,double inPhi
                          ,Particle *particleX
                          ,Particle *particleY
                         )

/* 

   For a two body decay, given cosTheta & Phi (GJ), & the total
   energy available (mass of the parent) this fills the
   four vectors of the daughter particles.

   Assumes we're in the rest frame of the parent.

*/

{

                              /* Make the 3-vectors we'll need */

Three b(gBeam->mFour.vec.space);
Three target(gTarget->mFour.vec.space);
Three recoil(gRecoil->mFour.vec.space);

Three n = recoil / target;
Three f = n / b;              /* slash = cross product if both
                                 arguments are vectors (they are) */

                              /* now unitize them */
b = b / b.Mag();
n = n / n.Mag();             /* just regular division */
f = f / f.Mag();

double SinTheta = 1 - SQR(inCosTheta);
       SinTheta = (SinTheta >= 0.0 ? sqrt(SinTheta) : 0.0 );

Three h    (  inCosTheta,
               SinTheta * cos(inPhi),
               SinTheta * sin(inPhi)  );
                

                              /* now do wacky math */
                              /* (Cramer's Rule)   */

SqMatrix3 A(b,n/b,f/b);

double detA = A.Determinant(); 

double detA1 = (A.ReplaceColumn(1,h)).Determinant(); 
double detA2 = (A.ReplaceColumn(2,h)).Determinant(); 
double detA3 = (A.ReplaceColumn(3,h)).Determinant(); 

Three answer(detA1/detA
            ,detA2/detA
            ,detA3/detA);  /* 'answer' holds lab frame angles */

/* now fill the two four vectors */

double S       = SQR(inTotalEnergy);
double energyX = (S + SQR(particleX->mMass) - SQR(particleY->mMass) ) /
                             (2.0 * inTotalEnergy);
double energyY  = inTotalEnergy - energyX;
double momentum = sqrt (SQR(energyX) - SQR(particleX->mMass));

answer = answer * momentum;  /* breaks total momentum into components */

/* finally, fill the four vectors */

particleX->mFour.Set( answer.vec.x, answer.vec.y, answer.vec.z,energyX);
particleY->mFour.Set(-answer.vec.x,-answer.vec.y,-answer.vec.z,energyY);
 
}

//////////////////////////////////////////////////////////////////

/*
 calculates CM energies for particleX and Y, fills them.
 You supply the total energy available.
*/


void Kinmaster::FillCMVectors(  
                           double inTotalEnergy
                          ,double inCosTheta
                          ,double inPhi
                          ,Particle *particleX
                          ,Particle *particleY
                         )
{

    /* find energies */

double S = SQR(inTotalEnergy);

double energyX = (S + SQR(particleX->mMass) - SQR(particleY->mMass) ) /
                             (2.0 * inTotalEnergy);
double energyY = inTotalEnergy - energyX;

    /* get momentum and components of momentum */

double momentum = sqrt (SQR(energyX) - SQR(particleX->mMass));

double pz = momentum * inCosTheta;

double SinTheta = 1.0 - SQR(inCosTheta);
      SinTheta = (SinTheta > 0.0 ? sqrt(SinTheta) : 0.0 );

double px = momentum * SinTheta * cos(inPhi);
double py = momentum * SinTheta * sin(inPhi);
    
    /* now fill the four vectors */

particleX->mFour.Set(px,py,pz,energyX);
particleY->mFour.Set(-px,-py,-pz,energyY);

}


//////////////////////////////////////////////////////////////////

/* 
   choose the momentum transfer variable "t" according to the user's
   specification.
*/

double Kinmaster::ChooseT( double energyProduct
                         ,double pmagProduct
                         ,Particle *beam
                         ,Particle *product)
{

  // we didn't use product->Momentum() or Energy() because the product four
  // vector has not been set (only the mass).  In fact, that's why we're 
  // finding t in this routine... 
  // so we can set the product and recoil four vectors.

  mTMin = ABS(SQR(beam->mMass) + SQR(product->mMass) - 
          2*(beam->Energy())*(energyProduct) + 
          2*(beam->Momentum())*(pmagProduct));

  mTMax = ABS(SQR(beam->mMass) + SQR(product->mMass) - 
          2*(beam->Energy())*(energyProduct) - 
          2*(beam->Momentum())*(pmagProduct));


  switch(mTransfer->transferCode)
    {
    case TDIST_1:
      mCurrentT = mTMin + (log(unfRand()) / (-(mTransfer->exponent))); 
      return (mCurrentT *= -1.0);
    case TDIST_2:
      while(TDistTwo((mCurrentT = Kinmaster::ChooseGuessT())) < unfRand());
      return (mCurrentT *= -1.0);
    case TDIST_3:
      return (mCurrentT = mUserTPicker->Pick()*(-1.0));
    case TDIST_USER:
      while((mCurrentT=mUserTPicker->Pick())<0.0);
      return (mCurrentT =(mTMin + mCurrentT)*(-1.0));
    case TDIST_DELTA:
      mCurrentT = mTransfer->t;
      return -(mCurrentT);   // user inputs abs(t)
    case TDIST_FLAT:
      mCurrentT = (mTransfer->t - mTransfer->tWidth) + 
                      ( unfRand() * ( 2.0 * mTransfer->tWidth) );
      return -(mCurrentT);
       
    default:
      fprintf(stderr, "Invalid T distribution specification\n");
      fprintf(stderr, "Check your input parameter file.\n\n");
      return 0;
    }
}

//////////////////////////////////////////////////////////////////

inline double Kinmaster::ChooseGuessT(void)
{
  // This function chooses an initial guess for t between tMin and tMax

  return mTMin + unfRand()*(mTMax-mTMin);
}


//////////////////////////////////////////////////////////////////

inline double Kinmaster::TDistOne(double t)
{
  return exp((double)(mTransfer->exponent)*t*(-1.0));
}

//////////////////////////////////////////////////////////////////

inline double Kinmaster::TDistTwo(double t)
{
  return E*(double)(mTransfer->exponent)*t*
         exp((double)(mTransfer->exponent)*t*(-1.0));
}

//////////////////////////////////////////////////////////////////

inline double Kinmaster::TDistThree(double t)
{
  return pow(t, (double)mTransfer->tPower)*
	 pow(t-mTMin, (double)mTransfer->spinPower)*
 	 exp((double)mTransfer->exponent*t*(-1.0))/
	 SQR(t+(double)SQR(mTransfer->exchangeMass));
}

//////////////////////////////////////////////////////////////////

inline double Kinmaster::BremThetaFunc(double theta)
{
	return theta / SQR(SQR(mThetaC) + SQR(theta));
}

//////////////////////////////////////////////////////////////////


void Kinmaster::FillTDistThree(double upperLimit, double binWidth)
{
  for(int i = 0; i < upperLimit/binWidth; ++i)
    bf1(TDIST3_HIST, binWidth*i+binWidth/2, TDistThree(binWidth*i+binWidth/2));
}

//////////////////////////////////////////////////////////////////
void Kinmaster::FillBremDist(double upperLimit, double binWidth)
{
  for(int i = 0; i < upperLimit/binWidth; ++i)
    bf1(BREMTHETA_HIST, binWidth*i+binWidth/2, BremThetaFunc(binWidth*i+binWidth/2));
}

//////////////////////////////////////////////////////////////////

double Kinmaster::GetProductCMEnergy(double TotCMEnergy,
                                    Particle *product,
                                    Particle *recoil)
{

double S = SQR(TotCMEnergy);
double productEnergy = (S + SQR(product->mMass) - SQR(recoil->mMass))/
                           (2*TotCMEnergy);
return productEnergy;

}

//////////////////////////////////////////////////////////////////


double Kinmaster::DoNumIntegration(void)  // Integrates t-dist 3
{
  double xLow = 0.0;
  const double stepSize = 0.01;
  double xHi = stepSize+xLow;
  double area = 0.0;


  double yHi, yLow, lastArea, convergence;

  do
    {
      yHi = (TDistThree)(xHi);
      yLow = (TDistThree)(xLow);
      
      lastArea = area;
      area += stepSize*MIN(yHi, yLow) + 0.5*stepSize*ABS(yHi-yLow);
           // square plus triangle

      convergence = ABS(lastArea-area)/lastArea;
      if(!finite(convergence))
	convergence = 50.0;

      xHi += stepSize;
      xLow += stepSize;
    }while(convergence > 0.000001);

  return xHi-stepSize;
}

//////////////////////////////////////////////////////////////////

void Kinmaster::LorentzXForm1(Four *cmFourVec, Four *source)
{

/* transform (source) to the frame where (cmFourVec) is at rest */

dvector3_t Beta;
double gamma = 0.0;
double mass= 0.0;

        Beta.x = cmFourVec->vec.space.x / cmFourVec->vec.t;
        Beta.y = cmFourVec->vec.space.y / cmFourVec->vec.t;
        Beta.z = cmFourVec->vec.space.z / cmFourVec->vec.t;

		mass = cmFourVec->mag4();

               // gamma equals E over m

        gamma = cmFourVec->vec.t / mass;
        

        double pdotbeta =
               (source->vec.space.x * Beta.x ) +
               (source->vec.space.y * Beta.y ) +
               (source->vec.space.z * Beta.z )  ;

        double factor =
               gamma*( (gamma*pdotbeta)/(gamma+1.0) -
                                      source->vec.t );

        Four change( Beta.x * factor
                    ,Beta.y * factor
                    ,Beta.z * factor
                    ,0.0
                   );

        double oldEnergy = source->vec.t;

        *source += change;
        source->vec.t = gamma*(oldEnergy - pdotbeta);


}
