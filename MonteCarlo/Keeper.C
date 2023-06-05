
/* test message */


// Keeper class code 

#include <fstream>
#include <iostream>
using namespace std;

#include <assert.h>
#include <stdio.h>
#include <math.h>

extern "C"{
#include <wabbit.h>
}
#include "Particle.h"
#include "util.h"
#include "targets.h"

#ifndef _KEEPER_CLASS_INCLUDED
 #include "Keeper.h"
#endif

#include "Kinmaster.h"


extern mc_param_t *gInParam;
extern particles_db_t *gDBPtr;
extern Kinmaster *gKinmaster;
extern Particle *gProduct;

#ifndef DO_E852
float lgd_dist = 710.8;
#endif

Keeper::Keeper(void)             // Keeper class constructor.  Initializes
{                              // number of particles to zero.
  mNumParts = 0;
}


void Keeper::keep(Particle *one)      // Gives Keeper the address of a Particle
{                                   // to be kept.  Checks to see if the max-
  assert(mNumParts < max_kept);     // imum number of Particles has been 
  mData[mNumParts] = one;           // exceded.
  ++mNumParts;
}


void Keeper::InitHistos(void)        // Sets up histograms 
{
//            Generated events
//
  bbook1
   (1+GENERATED_HISTOS,"Four photon eff mass GEN (GeV)", 100, 0.0, 3.0, 0);
  bbook1
   (2+GENERATED_HISTOS,"E individual photons GEN (GeV)", 100, 0.0, 20.0, 0);
  bbook1
   (3+GENERATED_HISTOS,"photon X (MPS-cm) at LGD GEN", 100, -152.0, 125.0, 0);
  bbook1
   (4+GENERATED_HISTOS,"photon Y (MPS-cm) at LGD GEN", 100, -90.0, 90.0, 0);
  bbook1
   (5+GENERATED_HISTOS,"t-distribution GEN", 251, -0.01, 2.5, 0.0);
  bbook1
   (6+GENERATED_HISTOS,"CosTh-GJ of Product decay GEN", 220, -1.1, 1.1, 0.0);
  bbook1
   (7+GENERATED_HISTOS,"Phi-GJ of Product decay GEN", 50, 0.0, 2*PI, 0.0);
  bbook1
   (8+GENERATED_HISTOS,"Daughter[0] mass (diagnostic) GEN", 250, 0.0, 2.5, 0.0);

//
//            Accepted events
//
  bbook1
   (1+ACCEPTED_HISTOS, "Four photon eff mass ACC (GeV)", 100, 0.0, 3.0, 0);
  bbook1
   (2+ACCEPTED_HISTOS, "E individual photons ACC (GeV)", 100, 0.0, 20.0, 0);
  bbook1
   (3+ACCEPTED_HISTOS, "photon X (MPS-cm) at LGD ACC", 100, -152.0, 125.0, 0);
  bbook1
   (4+ACCEPTED_HISTOS, "photon Y (MPS-cm) at LGD ACC", 100, -90.0, 90.0, 0);
  bbook1
   (5+ACCEPTED_HISTOS, "t-distribution ACC", 251, -0.01, 2.5, 0.0);
  bbook1
   (6+ACCEPTED_HISTOS, "CosTh-GJ of Product decay ACC", 100, 0.0, 1.0, 0.0);
  bbook1
   (7+ACCEPTED_HISTOS, "Phi-GJ of Prodcut decay ACC", 50, 0.0, 2*PI, 0.0);
  bbook1
   (8+ACCEPTED_HISTOS, "Daughter[0] mass (diagnostic) ACC", 250, 0.0, 2.5, 0.0);
  
}

////////////////////////////////////////////////////////////////////////////////

void Keeper::ConstructRotMatrix(double ctxz, double ctyz, double stxz, double styz)
{

	// this matrix takes vecs from the GJ frame to the overall CM frame
		
	Three row1(  ctxz,   -stxz*styz, ctyz*stxz  );
	Three row2(     0, 		   ctyz, styz       );
	Three row3( -stxz,   -ctxz*styz, ctxz*ctyz  );
	
	mRotMatrix.Set(row1,row2,row3);
	

}

////////////////////////////////////////////////////////////////////////////////

#ifdef DO_E852
unsigned int Keeper::e852cuts(void)      // Checks to see if an event passes
{                                      // the acceptance cuts.
  static position pho;
  static position tpho;
  int ret;

  static Box lgd(lgd_left, lgd_right, lgd_top, lgd_bottom);  
#ifdef DO_E852
  static Box dea(dea_left, dea_right, dea_top, dea_bottom);
#endif
  static Box beam(beam_left, beam_right, beam_top, beam_bottom);

  mAcceptBits = 0;  /* mAcceptBits holds the global cut results. */
                    /* 0 = good, 1 = bad */
                    /* See cutmasks.h */

  for (int i = 0; i < mNumParts; ++i)
    {
      mData[i]->mPersonalAcceptBits = 0;    /* reset individual accept bits */
      if (mCutCodes->doEnergySmear)         /* first smear the energy */
	SmearE(i);

      getPos(&pho, i); /* get position of particle at glass, csi, dea */

 /* SEPARATION CUTS */

      for (int j = i + 1; j < mNumParts; ++j)
	{
	  getPos(&tpho, j);
	  if (CircCut(gInParam->codes.separationCut
                     , pho.lgd.x
                     , pho.lgd.y
                     , &tpho.lgd))

            /* separation cut failed */
          {
           mData[i]->mPersonalAcceptBits =
              mData[i]->mPersonalAcceptBits | sepcut_mask;
           mData[j]->mPersonalAcceptBits =
              mData[j]->mPersonalAcceptBits | sepcut_mask;
          }; /* end if */
	}

 /* MINIMUM ENERGY CUT */

      if (mData[i]->mFour.vec.t < gInParam->codes.energyCut)
	mData[i]->mPersonalAcceptBits =
                 ( mData[i]->mPersonalAcceptBits | ecut_mask );

 /* LGD GEOMETRY CUT */

      if (lgd.RectCut(&(pho.lgd)))
	mData[i]->mPersonalAcceptBits =
                 ( mData[i]->mPersonalAcceptBits | geocut_mask );

 /* DEA GEOMETRY CUT */

      if ((ret=dea.RectCut(&(pho.dea))) != 0){
	mData[i]->mPersonalAcceptBits =
                 ( mData[i]->mPersonalAcceptBits | deacut_mask );
	switch(ret){
	case 1:
	  mData[i]->mPersonalAcceptBits =
	    ( mData[i]->mPersonalAcceptBits | deacut_mask_top );
	  break;
	case 2:
	  mData[i]->mPersonalAcceptBits =
	    ( mData[i]->mPersonalAcceptBits | deacut_mask_bottom );
	  break;
	case 3:
	  mData[i]->mPersonalAcceptBits =
	    ( mData[i]->mPersonalAcceptBits | deacut_mask_minusX );
	  break;
	case 4:
	  mData[i]->mPersonalAcceptBits =
	    ( mData[i]->mPersonalAcceptBits | deacut_mask_plusX );
	  break;
	default:
	  break;
	}
      }

 /* CSI GEOMETRY CUT */

      if (!CircCut(rcsi, e852_target_center_x, e852_target_center_y, &(pho.csi)))
	mData[i]->mPersonalAcceptBits =
                 ( mData[i]->mPersonalAcceptBits | csicut_mask );

 /* BEAM HOLE CUT */

      if (!beam.RectCut(&(pho.beam))) // beam hole;
	mData[i]->mPersonalAcceptBits =
                 ( mData[i]->mPersonalAcceptBits | beamcut_mask );

      /* now update the global accept bits */

      mAcceptBits = mAcceptBits | mData[i]->mPersonalAcceptBits;

      /* smear positions at LGD if necessary */

      if (mCutCodes->doPositionSmear)
	SmearPosition(&pho, i);
    }

  return mAcceptBits; /* this function returns the global AcceptBits */
}

#endif

////////////////////////////////////////////////////////

// plot the basic diagnostic histograms
// the argument histoBase is added on to the histogram index
//  for the purpose of plotting generated/accepted quantities

void Keeper::FillHistos(int histoBase)       // Fills histograms.
{

static position pho;
double efmass = 0.0;
static Four efMassVec;

  efMassVec.Set(0.0,0.0,0.0,0.0);
	
  for (int index = 0; index < mNumParts; ++index)
                    //index is the Particle number
  {
     // photon energies
 
     bf1(2+histoBase,(mData[index])->Energy(),1.0);

     // photon X,Y at the LGD

      getPos(&pho, index);

     bf1(3+histoBase,pho.lgd.x,1.0);
     bf1(4+histoBase,pho.lgd.y,1.0);

     efMassVec += (mData[index])->mFour;
  }; 

     // total photon effective mass 

    efmass = efMassVec.mag4();
 
    bf1(1+histoBase,efmass,1.0);

    // now fill the t distrution
    bf1(5+histoBase, (-1.0)*(gKinmaster->mCurrentT), 1.0);

    // fill costheta and phi for diagnostics
    bf1(6+histoBase, ABS(gProduct->mCosTheta), 1.0);
    bf1(7+histoBase, gProduct->mPhi, 1.0);

    // temporary file, should be removed later
    bf1(8+histoBase, gProduct->mDaughter[0]->mMass, 1.0);

}

////////////////////////////////////////////////////////

void Keeper::WriteHistos(void)      // Outputs histograms to a file.
{
  bbsave(0, "histos.bbook");
}


void Keeper::WriteHeader(char *filename)      // Opens output file, writes
{                                           // input parameters.
  FILE *temp = fopen(filename, "a");
  fwrite(gInParam, sizeof(*gInParam), 1, temp);
  fwrite(gDBPtr, sizeof(*gDBPtr), 1, temp);
}



////////////////////////////////////////////////////////

// Checks to see if point duo is within rad distance of point (orix, oriy)

int Keeper::CircCut(double rad, double orix, double oriy, point *duo)
{                                       
  return ((SQR(duo->x - orix) + SQR(duo->y - oriy)) < SQR(rad));
}


void Keeper::getPos(position *a, int index)     // Finds where particle hits
{                                             // various parts of the detector.
  double slope;                  //index = particle index

  slope = mData[index]->mFour.vec.space.x / mData[index]->mFour.vec.space.z;
  a->lgd.x = slope * (lgd_dist - mData[index]->mOrigin.z)+mData[index]->mOrigin.x;
#ifdef DO_E852
  a->dea.x = slope * (dea_dist - mData[index]->mOrigin.z)+mData[index]->mOrigin.x; 
  a->csi.x = slope * (csi_dist - mData[index]->mOrigin.z)+mData[index]->mOrigin.x;
#endif
  a->beam.x = a->lgd.x;
  slope = mData[index]->mFour.vec.space.y / mData[index]->mFour.vec.space.z;
  a->lgd.y = slope * (lgd_dist - mData[index]->mOrigin.z)+mData[index]->mOrigin.y;
#ifdef DO_E852
  a->dea.y = slope * (dea_dist - mData[index]->mOrigin.z)+mData[index]->mOrigin.y;
  a->csi.y = slope * (csi_dist - mData[index]->mOrigin.z)+mData[index]->mOrigin.y; 
#endif
  a->beam.y = a->lgd.y;
}


ostream &Keeper::print_data(ostream &aout)       // Prints out every Particle 
{                                              // in Keeper to stream aout.
  for(int i = 0; i < mNumParts; ++i)
    aout << *(mData[i]);
  return aout;
}

void Keeper::print_data(char *filename)        // Prints out every Particle in
{                                            // Keeper to file filenam.
  ofstream aout(filename,ios::app);
  for(int i = 0; i < mNumParts; ++i)
    aout << *mData[i];
}


void Keeper::cprint(FILE *cprnt, unsigned int accept)
{
  /* this routine write the info from the current event out to a file. */
  /* files are read and processed with separate analysis codes.        */

  c_file_header cheddar;  /* these types are defined in Keeper.h */
  c_file_data	bleu;
  
  vector3_t     vertex;   /* the priciple reaction vertex */

  vertex.x = gProduct->mOrigin.x;
  vertex.y = gProduct->mOrigin.y;
  vertex.z = gProduct->mOrigin.z;

  cheddar.numVectors = mNumParts;
  cheddar.globalAcceptBits = accept;
  cheddar.principleVertex = vertex;

  fwrite(&cheddar, sizeof(c_file_header), 1, cprnt);  /* write the header */

  for(int i = 0; i < mNumParts; ++i) {                /* write the data */
    bleu.vec = mData[i]->mFour.vec;
    bleu.PersonalAcceptBits = mData[i]->mPersonalAcceptBits;

    fwrite(&bleu, sizeof(c_file_data), 1, cprnt);
  }
}

/////////////////////////////////////////////////////////

// TransformAll
//
// we're passed a four-vector which defines the CM frame
// to transform to.  We transform all the particles into
// this frame.


void Keeper::TransformAll(Four *cmFourVec) 
{

dvector3_t Beta;
double gamma = 0.0;
double mass= 0.0;

         	// Beta equals p over E
         	// We need a three vector class!

	Beta.x = cmFourVec->vec.space.x / cmFourVec->vec.t;   
	Beta.y = cmFourVec->vec.space.y / cmFourVec->vec.t;
	Beta.z = cmFourVec->vec.space.z / cmFourVec->vec.t;

	mass = cmFourVec->mag4();

         	// gamma equals E over m

	gamma = cmFourVec->vec.t / mass;
        
         	// now loop over particles and transform them

	for(int index = 0; index < mNumParts; index++) 
                        // index is particle index
	{
	  double pdotbeta =
               ((mData[index])->mFour).vec.space.x * Beta.x + 
               ((mData[index])->mFour).vec.space.y * Beta.y + 
               ((mData[index])->mFour).vec.space.z * Beta.z   ; 

          double factor =
               gamma*( (gamma*pdotbeta)/(gamma+1) - 
                                      ((mData[index])->mFour).vec.t );

          Four change( Beta.x * factor
                      ,Beta.y * factor
		      ,Beta.z * factor
                      ,0.0
                     );
 
          (mData[index])->mFour += change;

          double oldEnergy = ((mData[index])->mFour).vec.t;
          ((mData[index])->mFour).vec.t = gamma*(oldEnergy - pdotbeta);
   

	}; /* done transforming particles */ 	
    
}

/////////////////////////////////////////////////////////

void Keeper::RotateGJtoCM(void)
{
   // Take all vecs from GJ frame to overall CM frame
   
   	for(int index = 0; index < mNumParts; index++) 
                        // index is particle index
	{
	 mData[index]->mFour = mRotMatrix * mData[index]->mFour;
	};
}

/////////////////////////////////////////////////////////

double Keeper::GetTotalEnergy(void)      // Returns sum of energy in every 
{                                      // Particle in Keeper.

double energy = 0.0;


	for(int index = 0; index < mNumParts; index++) 
                        // index is particle index
	{
	 energy += (mData[index])->Energy();
	}

	return energy;
}



// SmearE changes the energy of the particle; random with a gaussian
// distribution centered on the actual value.  Sigma is the resolution
// of the detector.


void Keeper::SmearE(int index)
{
  double xslope, yslope, zslope;
  xslope = mData[index]->mFour.vec.space.x / mData[index]->mFour.vec.t;
  yslope = mData[index]->mFour.vec.space.y / mData[index]->mFour.vec.t;
  zslope = mData[index]->mFour.vec.space.z / mData[index]->mFour.vec.t;

  double sigma = 
    mCutCodes->energySmearFloor + 
     mCutCodes->energySmearStat/(sqrt(mData[index]->mFour.vec.t));

  sigma *= mData[index]->mFour.vec.t;

  short done = 0;
  while (!done)   /* pick a new energy, making sure it's greater than 0 */
  {

    double newEnergy =  
     GetGaussianDist(mData[index]->mFour.vec.t, sigma);

    if (newEnergy > 0.0)
    {
     mData[index]->mFour.vec.t = newEnergy;
     done = 1;
    };

  } /* end while loop */

  mData[index]->mFour.vec.space.x = xslope * mData[index]->mFour.vec.t;
  mData[index]->mFour.vec.space.y = yslope * mData[index]->mFour.vec.t;
  mData[index]->mFour.vec.space.z = zslope * mData[index]->mFour.vec.t;  

//  cout << "E2 = " << SQR(mData[index]->mFour.vec.t) << " \n"; 
}



// SmearPosition changes the position of a particle to account for the 
// resolution of the detector.


void Keeper::SmearPosition(position *a, int index)
{
  double xslope, yslope, zslope, leng;
  a->lgd.x = GetGaussianDist(a->lgd.x, mCutCodes->positionSmear[0]);
  a->lgd.y = GetGaussianDist(a->lgd.y, mCutCodes->positionSmear[1]);
  leng = sqrt(SQR(mData[index]->mOrigin.x - a->lgd.x) + SQR(mData[index]->mOrigin.y - a->lgd.y) + SQR(lgd_dist - mData[index]->mOrigin.z));
  xslope = (a->lgd.x - mData[index]->mOrigin.x) / leng;
  yslope = (a->lgd.y - mData[index]->mOrigin.y) / leng;
  zslope = (lgd_dist - mData[index]->mOrigin.z) / leng;
  mData[index]->mFour.vec.space.x = xslope * mData[index]->mFour.vec.t;
  mData[index]->mFour.vec.space.y = yslope * mData[index]->mFour.vec.t;
  mData[index]->mFour.vec.space.z = zslope * mData[index]->mFour.vec.t;
}



//  mySmear same as SmearPosition.


void Keeper::mySmear(position *a, int index)
{
  double xslope, yslope;
  a->lgd.x = GetGaussianDist(a->lgd.x, mCutCodes->positionSmear[0]);
  a->lgd.y = GetGaussianDist(a->lgd.y, mCutCodes->positionSmear[1]);
  xslope = (a->lgd.x - mData[index]->mOrigin.x)/(lgd_dist - mData[index]->mOrigin.z);
  yslope = (a->lgd.y - mData[index]->mOrigin.y)/(lgd_dist - mData[index]->mOrigin.z);
  mData[index]->mFour.vec.space.z = mData[index]->mFour.vec.t / sqrt(SQR(xslope) + SQR(yslope) + 1);
  mData[index]->mFour.vec.space.x = xslope * mData[index]->mFour.vec.space.z;
  mData[index]->mFour.vec.space.y = yslope * mData[index]->mFour.vec.space.z;
}

