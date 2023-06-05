// Particle class code
//
// PCF 19 June 1995


#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef SOLARIS
#include <ieeefp.h>
#endif

#include "Keeper.h"
#include "Particle.h"
#include "Kinmaster.h"
#include "Three.h"
#include "util.h"

#define YLM_FILE "/ylms.bbook"
#define DO_HEL_SELECT 

static char *ylmFile;
extern Four into_rest_frame_of(Four vec1,Four vec2);
extern float angle_4v(Four vec1,Four vec2);

extern Four      *gLabFrameVec;   // allows us to know parent in lab frame
extern Keeper    *gMotherKeeper;  // she'll Lorentz transform for us
extern Particle  *gBeam;
extern Particle  *gTarget;        // need to define overall CM frame in decay2()
extern Particle  *gProduct;
extern Particle  *gRecoil;
extern Kinmaster *gKinmaster;     // everyone loves the Kinmaster!!
                                  // takes care of filling in daughter
                                  // particle four vectors.

////////////////////////////////////////////////////

Particle::Particle() // constructor
{
  mOrigin.x = 0.0;
  mOrigin.y = 0.0;
  mOrigin.z = 0.0;

  for(int j = 0; j < max_daughters; ++j)
    mDaughter[j] = NULL;

  mMass = 0.0;
  mFour = Four();
  mNumDaughters = 0;
  mParent = NULL;


  /* allocate space for the NosPickers */
  
  mUserMassPicker = NULL;
  mUserMassPicker = new(NosPicker);
  assert(mUserMassPicker != NULL);

  mUserCosThetaPicker = NULL;
  mUserCosThetaPicker = new(NosPicker);
  assert(mUserCosThetaPicker != NULL);

  mUserPhiPicker = NULL;
  mUserPhiPicker = new(NosPicker);
  assert(mUserPhiPicker != NULL);

  mI=NULL;
  mScale=NULL;
}

Particle::~Particle() // Destructor
{
}

////////////////////////////////////////////////////

void Particle::InitParticle()
{
  /* This routine sets up the machinery that picks random numbers
     from a user specified distribution.
     The user may be supplying cosTheta and phi distributions in 
     a histogram.
     Also, for Y_LM decays, we have histograms already on file that
     need to be read.
  */
 
  /* init the NosPickers

     There is a histogram file which contains
     all of the ylms that we are using. 
     Run the "createYLM" program to create this file.
     The histoID for each histogram in this file 
     (ylms.bbook) is stored as a symbolic constant in Particle.h
  */

  if(!ylmFile){
    char *dir;
    char *subdir="/MonteCarlo/";

    dir = getenv("RADPHI_CONFIG");
    assert(dir);
    ylmFile = (char *)calloc(strlen(dir)+strlen(subdir)+strlen(YLM_FILE)+1,sizeof(char));
    strcpy(ylmFile,dir);
    strcat(ylmFile,subdir);
    strcat(ylmFile,YLM_FILE);
  }	    

  if(mInfo->massCode == MASS_USER)
  {
    if(mUserMassPicker->InitNosPicker(mInfo->massHisto.id, 
				      mInfo->massHisto.filename))
      {
	/* error message handling has already been done */
	fprintf(stderr, "Check your input parameter file for for possible");
	fprintf(stderr, "sources of error\n");
	fprintf(stderr, "Exiting...\n\n");
	exit(EXIT_FAILURE);
      };
   }; /* end mass init */

  
  if(((mInfo->twoBody.code & DECAY_THETA_MASK) == DECAY_THETA_CORRELATED)){
    
    /* read in a 3-d (m,costheta,phi) intensity distribution */
    
    if(ReadIntensityFile(mInfo->twoBody.thetaHisto.filename,&mI,
			 &mNMassBins, &mMassLow, &mMassBinSize, &mNCthBins, &mNPhiBins, mScale))
      {
	fprintf(stderr,"\nfailed to read intensity file ... bailing out");
	exit(EXIT_FAILURE);    
      };
    
    
  }else{
    
    if((mInfo->twoBody.code & DECAY_THETA_MASK) == DECAY_THETA_USER)
      {
	if(mUserCosThetaPicker->InitNosPicker(mInfo->twoBody.thetaHisto.id, 
					      mInfo->twoBody.thetaHisto.filename))
	  {
	    /* error message handling has already been done */
	    fprintf(stderr, "Check your input parameter file for for possible");
	    fprintf(stderr, "sources of error\n");
	    fprintf(stderr, "Exiting...\n\n");
	    exit(EXIT_FAILURE);
	  };
      }; /* end theta init */
    
    
    if((mInfo->twoBody.code & DECAY_PHI_MASK) == DECAY_PHI_USER)
      {
	if(mUserPhiPicker->InitNosPicker(mInfo->twoBody.phiHisto.id,
					 mInfo->twoBody.phiHisto.filename))
	  {
	    /* error message handling has already been done */
	    fprintf(stderr, "Check your input parameter file for for possible");
	    fprintf(stderr, "sources of error\n");
	    fprintf(stderr, "Exiting...\n\n");
	    exit(EXIT_FAILURE);
	  };
      }; /* end phi init */
  }
  
  if((mInfo->twoBody.code & DECAY_THETA_MASK) == DECAY_THETA_YLM)
  {
    if(mInfo->decay[0].nProducts)
      switch(mInfo->spin)
	{
	case 0:
	  break;
	  
	case 1:
	  if(mUserCosThetaPicker->InitNosPicker(HIST_Y10, ylmFile))
	    {
	      fprintf(stderr, "Problems finding histogram 10 in %s\n",ylmFile);
	      fprintf(stderr, "Run createYLM and try again\n");
	      fprintf(stderr, "Exiting...\n\n");
	      exit(EXIT_FAILURE);
	    }
	  break;
	  
	case 2:
	  if(mUserCosThetaPicker->InitNosPicker(HIST_Y20, ylmFile))
	    {
	      fprintf(stderr, "Problems finding histogram 20 in %s\n",ylmFile);
	      fprintf(stderr, "Run createYLM and try again\n");
	      fprintf(stderr, "Exiting...\n\n");
	      exit(EXIT_FAILURE);
	    }
	  break;
	  
      case 4:
	  if(mUserCosThetaPicker->InitNosPicker(HIST_Y40, ylmFile))
	    {
	      fprintf(stderr, "Problems finding histogram 40 in %s\n",ylmFile);
	    fprintf(stderr, "Run createYLM and try again\n");
	      fprintf(stderr, "Exiting...\n\n");
	    exit(EXIT_FAILURE);
	    }
	  break;

	default:
	  fprintf(stderr, "Sorry, angular decays of spin %d ", mInfo->spin);
	  fprintf(stderr, "are not currently implemented\n");
	  fprintf(stderr, "Exiting...\n\n");
	  exit(EXIT_FAILURE);
	}	
    }; /* end YLM init */
}


////////////////////////////////////////////////////

double Particle::Energy(void)
{
 return mFour.vec.t;
}

////////////////////////////////////////////////////

double Particle::Momentum(void)
{
 return mFour.mag3();
}

////////////////////////////////////////////////////

double Particle::Mass(void)
{
 return mFour.mag4();  /* just a hook to the four vector function */
}

////////////////////////////////////////////////////

double Particle::MassSq(void)
{
 return mFour.mag4Sq();  /* just a hook to the four vector function */
}

////////////////////////////////////////////////////

void Particle::Decay(  Keeper *userKeeper
                     , Keeper *photonKeeper
                     , int doGJ )
{
  if (mInfo->nDecays) 
  {
    
    switch (mNumDaughters)
    {
      case 2: decay2(doGJ); //decay2 sets up two daugther particles
   	      break;

      case 3: decay3(); //decay3 sets up three daughter particles
              break;

      default: // OOPS! we can't do it.
             printf("\n\n Sorry, one of the decays you're asking \n");
             printf("for is not currently implemented (Particle::Decay).\n\n");
             exit(EXIT_FAILURE);
    };
    
    for(int i = 0; i < mNumDaughters; ++i)
      (*mDaughter[i]).Decay(userKeeper,photonKeeper,0);  /* 0 means no GJ
                                                            inverse */ 

  };                                      //daughters decay


}

////////////////////////////////////////////////////

void Particle::SetEnergy(double newEnergy)
{
  mFour.vec.t = newEnergy;
}

////////////////////////////////////////////////////

void Particle::operator=(Particle &one)
{
  
  mInfo = one.mInfo;

  mOrigin.x = one.mOrigin.x;
  mOrigin.y = one.mOrigin.y;
  mOrigin.z = one.mOrigin.z;

  for(int j = 0; j < max_daughters; ++j)
    mDaughter[j] = one.mDaughter[j];
  mMass = one.mMass;
  mFour = one.mFour;

  mNumDaughters = one.mNumDaughters;
  mParent = one.mParent;
}

////////////////////////////////////////////////////



ostream &operator<<(ostream &aout, Particle &one)
{
  aout << one.mInfo->name << '\n';
  aout << "PX = " << one.mFour.vec.space.x << '\n';
  aout << "PY = " << one.mFour.vec.space.y << '\n';
  aout << "PZ = " << one.mFour.vec.space.z << '\n';
  aout << "E = " << one.mFour.vec.t << '\n' << '\n';
  return aout;
}


////////////////////////////////////////////////////

/* decay into two daugter particles */

void Particle::decay2(int doGJ)
{
  static Four overallCMVec;
  static Four tempBeam;
  static Four tempTarget;
  Four overallCM(0.0,0.0,0.0,1.0);
  int must_select_phi=1;

/*
// Particles start in the GJ frame.
// First transform ALL particles into the center of mass of this particle.
*/

  /*fprintf(stdout,"recoil CM before:\n");
  gRecoil->mFour.print(); */
  /*  Four gproddummy=into_rest_frame_of(gProduct->mFour,this->mFour);
  Four grecdummy=into_rest_frame_of(gRecoil->mFour,this->mFour);
  Four newCM=into_rest_frame_of(overallCM,this->mFour); */

  gKinmaster->LorentzXForm1(&(this->mFour),&overallCM);
  gMotherKeeper->TransformAll( &(this->mFour) ); 
  

  /*  fprintf(stdout,"product recoil CM after real transform:\n");
      gProduct->mFour.print();
  gRecoil->mFour.print(); 
  overallCM.print();
  fprintf(stdout,"product recoil CM after my transform:\n");
  gproddummy.print();
  grecdummy.print();
  newCM.print(); */



//	tempBeam = gBeam->mFour;
//	tempTarget = gTarget->mFour;  
	
	// put tempBeam and tempTarget into this particle's cm frame
	
//	gKinmaster->LorentzXForm1(&(this->mFour),&tempBeam);
//	gKinmaster->LorentzXForm1(&(this->mFour),&tempTarget);

/*
// Now pick cosTheta and phi, the angle at which the first daughter particle
// comes off in the parent's CM frame.
*/

  

  if(((mInfo->twoBody.code & DECAY_THETA_MASK) == DECAY_THETA_CORRELATED)){
     
    /* use this mFour to find the mass, then figure out which mbin 
     and the intensity function from the array to find a costheta and phi */
 
    int mbin=(int)((mMass-mMassLow)/mMassBinSize);

    /* some bounds check to avoid SIGSEV */
    if(mMass<mMassLow){
      mbin=0;
    }
    if(mbin>=mNMassBins){
      mbin=mNMassBins-1;
    }
    
    int cthBin=(int)((float)mNCthBins*unfRand()); 
    /* these will cast such that the fractional part is thrown away */
    int phiBin=(int)((float)mNPhiBins*unfRand()); 
    
    while(*(mI+mbin*mNCthBins*mNPhiBins+cthBin*mNPhiBins+phiBin)<unfRand()){
      /* I is normalized s.t. max value is 1.0 */
      cthBin=(int)((float)mNCthBins*unfRand());  
      phiBin=(int)((float)mNPhiBins*unfRand()); 
    }
    /* now assign the double precision value to costheta and phi */
    mCosTheta=-1.0 + ((float)cthBin+unfRand())*2.0/(float)mNCthBins;
    mPhi= ((float)phiBin+unfRand())*2.0*PI/(float)mNPhiBins;
    
  }
  else{
    
    /* first choose theta */
 
    /*    fprintf(stdout,"[debug:] twobody.code = %x.\n",mInfo->twoBody.code);
    exit(1); */
    switch(mInfo->twoBody.code & DECAY_THETA_MASK) 
      {
      case DECAY_THETA_UNIFORM:              /* isotropic decay */
	
	mCosTheta = -1.0 + unfRand()*2.0;
	break;
	
      case DECAY_THETA_FLAT:                 /* flat in user specf'd region */
	
	mCosTheta = mInfo->twoBody.theta - mInfo->twoBody.thetaWidth +
	  unfRand()*( 2* mInfo->twoBody.thetaWidth );
	break;
	
      case DECAY_THETA_YLM:
	switch(mInfo->spin)
	  {
	  case 0:                                     // same as flat decay 
	    mCosTheta = -1.0 + unfRand()*2.0;
	    break;
	  case 1:                                     // then they want a y10
	    mCosTheta = mUserCosThetaPicker->Pick();
	    break;
	  case 2:                                     // then they want a y20
	    mCosTheta = mUserCosThetaPicker->Pick();
	    break;
	  case 4:                                     // then they want a y40
	    mCosTheta = mUserCosThetaPicker->Pick();
	    break;
	  default:
	    fprintf(stderr, "Sorry, angular decays of spin %d ", mInfo->spin);
	    fprintf(stderr, "are not currently implemented\n");
	    fprintf(stderr, "Exiting...\n\n");
	    exit(EXIT_FAILURE);
	  }; /* end inner switch on spin */
	break;

      case DECAY_THETA_DALITZ:
	printf("\n\n You are asking for a two-body dalitz decay! \n");
	printf("Exiting...\n\n");
	exit(EXIT_FAILURE);
	break;
	
      case DECAY_THETA_REAL_PHI:
	mCosTheta=mUserCosThetaPicker->Pick_real_phi(gRecoil->mFour,
						     overallCM,
						     &mPhi);
	must_select_phi=0;
	break;

      case DECAY_THETA_USER:
	mCosTheta = mUserCosThetaPicker->Pick();
	break;
	
      default:
	fprintf(stderr, "Invalid angular decay specification.\n");
	fprintf(stderr, "Check your input parameter file.\n\n");
	exit(EXIT_FAILURE);
	
	break;
	
      }; /* end switch for finding cosTheta*/
    
    /* now choose phi (in radians) */
    
    if(must_select_phi){
      switch(mInfo->twoBody.code & DECAY_PHI_MASK) 
	{
	case DECAY_PHI_FLAT:  /* pick uniformly from user specf'd region */
	  mPhi = mInfo->twoBody.phi - mInfo->twoBody.phiWidth +
	    unfRand()*( 2* mInfo->twoBody.phiWidth );
	  
	  break;
	  
	case DECAY_PHI_UNIFORM:
	  mPhi = 2.0*PI*unfRand();
	  break;
	  
	case DECAY_PHI_USER: 
	  mPhi = mUserPhiPicker->Pick();
	  break;
	  
	};
      /* end phi switch */
    }
  }
  
// The total energy available is the mass of the parent particle.
// Pass this and the angular info along to the Kinmaster, who will
// fill in the daughter particles' four vectors.  The daughter particles
// masses are chosen when the parent's mass is chosen. 

   	double energyAvail = this->mMass;
//        double energyAvail = this->mFour.vec.t;     

             // should be the same as the mass--
             // different because of roundoff errors

 switch (doGJ)
 {
  case 1:  /* i.e., the main reaction product */

   gKinmaster->
     Fill2BodyGJInverse(energyAvail
                       ,this->mCosTheta
                       ,this->mPhi
                       ,this->mDaughter[0]
                       ,this->mDaughter[1]
                       );
                       
 		       break;

  default:
   gKinmaster->
     FillCMVectors(energyAvail
                       ,this->mCosTheta
                       ,this->mPhi
                       ,this->mDaughter[0]
                       ,this->mDaughter[1]
                       );
  
                       break;
                       
  }; /* end switch */


 

// Set the daugther particles' origin.  It's the same as the parent
// if the parent isn't stable.  Otherwise we propagate the parent's
// origin.

   static double sovertau;
   static double clength;
   static double magP;

   if (mInfo->lifetime != 0.0)  /* stable particle */
    {
     /* propagate the particle's origin to find daughters' origin */
    
    do
    {
     sovertau = -log(1.0 - unfRand());
    } while ( !(finite(sovertau)) );   /* make sure we don't NAN to hell */

    magP = mFour.mag3();  /* magnitude of our momentum */
      // BUG:
      // mFour.mag3() = 0 since we're still in the CM frame.
      // must x-form to lab frame first
    
    printf("\nWARNING: finite particle lifetime code not complete");
    printf("\nYour results will not be valid.\n\n");
      
    clength = (magP / mMass) * (mInfo->lifetime) * (sovertau);

    mOrigin.x += (mFour.vec.space.x / magP) * clength;
    mOrigin.y += (mFour.vec.space.y / magP) * clength;
    mOrigin.z += (mFour.vec.space.z / magP) * clength;

   }

   mDaughter[0]->mOrigin.x = mOrigin.x;
   mDaughter[1]->mOrigin.x = mOrigin.x;

   mDaughter[0]->mOrigin.y = mOrigin.y;
   mDaughter[1]->mOrigin.y = mOrigin.y;

   mDaughter[0]->mOrigin.z = mOrigin.z;
   mDaughter[1]->mOrigin.z = mOrigin.z;
  
// Now transform daughter particles back into the GJ frame.

  overallCMVec = gBeam->mFour + gTarget->mFour;
//    overallCMVec = tempBeam + tempTarget;
   gMotherKeeper->TransformAll(&overallCMVec);

//	gKinmaster->LorentzXForm1(&overallCMVec,&(mDaughter[0]->mFour));
//	gKinmaster->LorentzXForm1(&overallCMVec,&(mDaughter[1]->mFour));



}

////////////////////////////////////////////////////
void Particle::decay3()   /* spin 0 three body decay */
{
  static Four overallCMVec;
  
  double m, m1, m2, m3, m12, m13, m23;
  double e1, e2, e3, q1, q3, p1, p2, p3;
  double qf, min12, max12, delta12, delta23, 
         min23, max23, max13, min13, delta13;
  double cosa, cosb, sina, sinb, temp1, temp2;
  
  // Transform all particles into the decaying particle's center of mass frame.
  
  gMotherKeeper->TransformAll( &(this->mFour) );
  
  m = mMass;
  m1 = mDaughter[0]->mMass;
  m2 = mDaughter[1]->mMass;
  m3 = mDaughter[2]->mMass;
  
  qf = SQR(m1) + SQR(m2) + SQR(m3) + SQR(m);
  
  // Compute rough (square) Dalitz boundaries
  
  min12 = SQR(m1 + m2);
  max12 = SQR(m - m3);
  delta12 = max12 - min12;
  
  min23 = SQR(m2 + m3);
  max23 = SQR(m - m1);
  delta23 = max23 - min23;
  
  min13 = SQR(m1 + m3);
  max13 = SQR(m - m2);
  delta13 = max13 - min13;
  	
  do 
  {
    do 
    {
      do 
      {
		
		/* pick approximate values for m12, m23, m13 */
		  
	  		m12 = min12 + delta12 * unfRand();
	  		m23 = min23 + delta23 * unfRand();
	  		m13 = qf - m12 - m23;
		 		
  		
		/* Now do a detailed check to make sure m12, m23 are
		   inside the Dalitz boundaries.  */
		
		/* see 1994 Particle Data Book, p. 1291 */
		
		e1 = (m12 + SQR(m1) - SQR(m2)) / (2 * sqrt(m12));
		e3 = (SQR(m) - m12 - SQR(m3)) / (2 * sqrt(m12));
		q1 = sqrt(SQR(e1) - SQR(m1));
		q3 = sqrt(SQR(e3) - SQR(m3));
		
		max13 = SQR(e1 + e3) - SQR(q1 - q3);
		min13 = SQR(e1 + e3) - SQR(q1 + q3);
      
      } while ((m13 > max13) || (m13 < min13));
      
 
      /* Have masses; now compute energies and momenta */
      
      e3 = (SQR(m) + SQR(m3) - m12) / (2 * m);
      e2 = (SQR(m) + SQR(m2) - m13) / (2 * m);
      e1 = m - e2 - e3;
      
      p1 = sqrt(SQR(e1) - SQR(m1));
      p2 = sqrt(SQR(e2) - SQR(m2));
      p3 = sqrt(SQR(e3) - SQR(m3));
      
      cosb = (SQR(p1) + SQR(p3) - SQR(p2)) / (2 * p1 * p3);
    
    } while ((cosb > 1.0) || (cosb < -1.0));
    
    sinb = sqrt(1.0 - SQR(cosb));
    cosa = (p1 - p3 * cosb) / p2;
  
  } while ((cosa > 1.0) || (cosa < -1.0));
  
  
  /* Now randomize the orientation of the event. */
  
  sina = sqrt(1.0 - SQR(cosa));
  double phi = 2.0 * PI * unfRand();
  double cp = cos(phi);
  double sp = sin(phi);
  double ct = 2 * unfRand() - 1.0;
  double st = sqrt(1.0 - SQR(ct));
  double alf = 2.0 * PI * unfRand();
  double ca = cos(alf);
  double sa = sin(alf);
  
  /* Compute rotation matrix elements */
  
  double r11 = (ca * ct * cp) - (sa * sp);
  double r13 = ca * st;
  double r21 = - (sa * ct * cp) - (ca * sp);
  double r23 = - sa * st;
  double r31 = - st * cp;
  double r33 = ct;
  
  /* Finally, set the daughter particle four vectors. */
  
  mDaughter[0]->mFour.vec.space.x = r11 * p1;
  mDaughter[0]->mFour.vec.space.y = r21 * p1;
  mDaughter[0]->mFour.vec.space.z = r31 * p1;
  mDaughter[0]->mFour.vec.t = e1;
  
  temp1 = - p2 * cosa;
  temp2 = - p2 * sina;
  mDaughter[1]->mFour.vec.space.x = r11 * temp1 + r13 * temp2;
  mDaughter[1]->mFour.vec.space.y = r21 * temp1 + r23 * temp2;
  mDaughter[1]->mFour.vec.space.z = r31 * temp1 + r33 * temp2;
  mDaughter[1]->mFour.vec.t = e2;

  temp1 = - p3 * cosb;
  temp2 = p3 * sinb;
  mDaughter[2]->mFour.vec.space.x = r11 * temp1 + r13 * temp2;
  mDaughter[2]->mFour.vec.space.y = r21 * temp1 + r23 * temp2;
  mDaughter[2]->mFour.vec.space.z = r31 * temp1 + r33 * temp2;
  mDaughter[2]->mFour.vec.t = e3;

// Set the daugther particles' origin.  It's the same as the parent
// if the parent isn't stable.  Otherwise we propagate the parent's
// origin.

  static double sovertau;
  static double clength;
  static double magP;

  if (mInfo->lifetime != 0.0)  /* stable particle */
    {
      /* propagate the particle's origin to find daughters' origin */
      
      do
	{
	  sovertau = -log(1.0 - unfRand());
	} while ( !(finite(sovertau)) );   /* make sure we don't NAN to hell */
      
      magP = mFour.mag3();  /* magnitude of our momentum */
      
      clength = (magP / mMass) * (mInfo->lifetime) * (sovertau);
      
      mOrigin.x += (mFour.vec.space.x / magP) * clength;
      mOrigin.y += (mFour.vec.space.y / magP) * clength;
      mOrigin.z += (mFour.vec.space.z / magP) * clength;
      
    }

  for(int i = 0; i < 3; ++i) 
    {
      mDaughter[i]->mOrigin.x = mOrigin.x;
      mDaughter[i]->mOrigin.y = mOrigin.y;
      mDaughter[i]->mOrigin.z = mOrigin.z;
    }

// Return to overall CM frame

  overallCMVec = gTarget->mFour + gBeam->mFour;
  gMotherKeeper->TransformAll(&overallCMVec);

}

////////////////////////////////////////////////////
  
void Particle::ChooseMass(void)
{

  double testMass;      // randomly produced test mass used in BW and Gaussian 
  double massDaughters; // sum of the masses of the daughter particles
  int i;

  // This function searches the massCode in mInfo to determine what type of 
  // distribution to choose the mass from, then it stores that mass in mMass.
  // It also sets the particles 4-vector to the CM frame
  
  switch(mInfo->massCode)
    {
    case MASS_FLAT:
      do
	{
	  /* first set the mass of the daughter particles to 0 */
	  massDaughters = 0.0;
	  
	  /* now choose the mass of the daughters and sum the mass */
	  for(i = 0; i < mNumDaughters; ++i)
	    {
	      mDaughter[i]->ChooseMass();
	      massDaughters += mDaughter[i]->mMass;
	    }

	  /* now get the mass of the actual particle */
	  mMass = (mInfo->mass - (mInfo->sigma)) + unfRand()*(2*(mInfo->sigma));
	}
      while(mMass < massDaughters);

      mFour.Set(0.0, 0.0, 0.0, mMass);
      break;
    case MASS_BW:
      do
	{
	  /* first set the mass of the daughter particles to 0 */
	  massDaughters = 0.0;

          /* now choose the mass of the daughters and sum the mass */
	  for(i = 0; i < mNumDaughters; ++i)
	    {
	      mDaughter[i]->ChooseMass();
	      massDaughters += mDaughter[i]->mMass;
	    }
          
	  /* now get the mass of the actual particle */
	  while(breitWigner((double)mInfo->sigma, (double)mInfo->mass, 
	       (testMass = Particle::ChooseGuessMass())) < unfRand());

	  mMass = testMass;
	} 
      while(mMass < massDaughters);
 
      mFour.Set(0.0, 0.0, 0.0, mMass);
      break;
    case MASS_GAUSSIAN:
      do
	{
	 /* first set the mass of the daughter particles to 0 */
	  massDaughters = 0.0;

          /* now choose the mass of the daughters and sum the mass */
	  for(i = 0; i < mNumDaughters; ++i)
	    {
	      mDaughter[i]->ChooseMass();
	      massDaughters += mDaughter[i]->mMass;
	    }
          
	  /* now get the mass of the actual particle */
	  while(gaussian((double)mInfo->sigma, (double)mInfo->mass, 
	       (testMass = Particle::ChooseGuessMass())) < unfRand());

	  mMass = testMass;
	} 
      while(mMass < massDaughters);

      mFour.Set(0.0, 0.0, 0.0, mMass);
      break;
    case MASS_DELTA:
      do
	{
	  /* first set the mass of the daughter particles to 0 */
	  massDaughters = 0.0;

          /* now choose the mass of the daughters and sum the mass */
	  for(i = 0; i < mNumDaughters; ++i)
	    {
	      mDaughter[i]->ChooseMass();
	      massDaughters += mDaughter[i]->mMass;
	    }
          
	  /* now get the mass of the actual particle */
	  mMass = mInfo->mass;
	}
      while(mMass < massDaughters);

      mFour.Set(0.0, 0.0, 0.0, mMass);
      break;
    case MASS_USER:
      do
	{
	  /* first set the mass of the daughter particles to 0 */
	  massDaughters = 0.0;
	  
          /* now choose the mass of the daughters and sum the mass */
	  for(i = 0; i < mNumDaughters; ++i)
	    {
	      mDaughter[i]->ChooseMass();
	      massDaughters += mDaughter[i]->mMass;
	    }
          
	  /* now get the mass of the actual particle */
	  mMass = mUserMassPicker->Pick();
	}
      while(mMass < massDaughters);

      mFour.Set(0.0, 0.0, 0.0, mMass);
      break;
    default:
      fprintf(stderr, "Invalid mass distribution specification\n");
      fprintf(stderr, "Check your input parameter file\n\n");
      return;
    }
}

////////////////////////////////////////////////////

inline double Particle::ChooseGuessMass(void)
{

  // choose the mass of the particle with a low mass greater than the sum of
  // its daughters and a high mass within 4 sigma of the mean

  return ((double)(mInfo->mass - 3.0*mInfo->sigma) + 
	  unfRand()*((double)(6.0*mInfo->sigma)));
}
   

////////////////////////////////////////////////////

/* Four into_rest_frame_of(Four vec1,Four vec2)
{
  Four beta,transformed;
  double energy,betavalsq,bX,bY,bZ,gamma,gm_one;
  //  double temp=vec2.vec.t;
  beta=vec2*(1/vec2.vec.t);
  bX=beta.vec.space.x;
  bY=beta.vec.space.y;
  bZ=beta.vec.space.z;
  betavalsq=SQR(beta.mag3());
  gamma=vec2.vec.t*(1/vec2.mag4());
  gm_one=gamma-1.0E+0;
  transformed.vec.t=
    (gamma*vec1.vec.t)+
    (-(gamma*bX)*vec1.vec.space.x)+
    (-(gamma*bY)*vec1.vec.space.y)+
    (-(gamma*bZ)*vec1.vec.space.z);
  transformed.vec.space.x=
    (-(gamma*bX)*vec1.vec.t)+
    ((1+(gm_one*(SQR(bX)/betavalsq)))*vec1.vec.space.x)+
    (((gm_one*bX*bY)/betavalsq)*vec1.vec.space.y)+
    (((gm_one*bX*bZ)/betavalsq)*vec1.vec.space.z);
  transformed.vec.space.y=
    (-(gamma*bY)*vec1.vec.t)+
    (((gm_one*bX*bY)/betavalsq)*vec1.vec.space.x)+
    ((1+(gm_one*(SQR(bY)/betavalsq)))*vec1.vec.space.y)+
    (((gm_one*bY*bZ)/betavalsq)*vec1.vec.space.z);
  transformed.vec.space.z=
    (-(gamma*bZ)*vec1.vec.t)+
    (((gm_one*bX*bZ)/betavalsq)*vec1.vec.space.x)+
    (((gm_one*bY*bZ)/betavalsq)*vec1.vec.space.y)+
    ((1+(gm_one*(SQR(bZ)/betavalsq)))*vec1.vec.space.z);  
  return transformed;
  
}*/

/* float angle_4v(Four vec1,Four vec2)
{
  double one_dot_two,angle,arg;
  one_dot_two=(vec1.vec.space.x*vec2.vec.space.x)+(vec1.vec.space.y*vec2.vec.space.y)+(vec1.vec.space.z*vec2.vec.space.z);
  arg = one_dot_two/(vec1.mag3()*vec2.mag3());
  if(fabs(arg)>1){
    if (arg<0)
      return(PI);
    else
      return(0);
  }
  angle=acos(one_dot_two/(vec1.mag3()*vec2.mag3()));
  return angle;
} */
