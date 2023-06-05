// Input.C
//
// reads input file and sets up parameters for NewBnlgams
//
// PCF 21 June 1995
//

#include <assert.h>

#include "Input.h"
#include "Kinmaster.h"
#include "Keeper.h"
#include "Particle.h"



// links to various and sundry globals

extern mc_param_t 	*gInParam;
extern particles_db_t  	*gDBPtr;
extern Keeper 		*gMotherKeeper; 
extern Keeper 		*gUserKeeper; 
extern Keeper 		*gPhotonKeeper; 
extern Kinmaster 	*gKinmaster;
extern Particle 	*gBeam;
extern Particle 	*gTarget;
extern Particle 	*gProduct;
extern Particle 	*gRecoil;

/////////////////////////////////////////////////////////////////

Input::Input()
{

}

/////////////////////////////////////////////////////////////////

int Input::ReadInput(char *filename)
{
const int error = 1;

     	if(initParam(filename)) /* non zero if it fails */
         return error;
 int 	isReadOK = loadParam(gInParam,&gDBPtr);

 return isReadOK;
}

/////////////////////////////////////////////////////////////////

void Input::SetUpParticles(void)
{
// at this point, CreateGlobalObject has been called, but we still need
// to set up every particle beside the beam, target, product, and recoil.
// We also need to give pointers to the input information to various objects.

// first, give the beam, target, product, and recoil their database info
// (of particle_db_t defined in mc_param.h)

   gBeam->mInfo    = & ( gDBPtr->particle[gInParam->reaction.beam] );
   gTarget->mInfo  = & ( gDBPtr->particle[gInParam->reaction.target]  );

// now we have to call InitParticle, because some of the initialization of the 
// particle can't be carried out intil mInfo is set.  Basically this just 
// inits the NosPickers for cosTheta and mass in the case when the user has 
// as user specified distribution in the form of a histo

   gBeam->InitParticle();
   gTarget->InitParticle();

// Now give the Kinmaster the beam and momentum transfer info   

   gKinmaster->mBeam     = &( gInParam->beam );
   gKinmaster->mTransfer = &( gInParam->transfer );

// as with the particles, we now have to call InitKinmaster which does some 
// extra initialization after mTransfer is set.

   gKinmaster->InitKinmaster();

// Tell the keepers about the cut info

   gMotherKeeper->mCutCodes = &( gInParam->codes );
   gUserKeeper->mCutCodes   = &( gInParam->codes );
   gPhotonKeeper->mCutCodes = &( gInParam->codes );

// Assign the beam and target to the proper Keepers.
// The rest of the particles are handled in FillTree.

   gMotherKeeper->keep(gBeam);
   gMotherKeeper->keep(gTarget);  // MotherKeeper gets them all

    // UserKeeper gets them if the user wants them

   if (gBeam->mInfo->returnData) gUserKeeper->keep(gBeam);   
   if (gTarget->mInfo->returnData) gUserKeeper->keep(gTarget);

// Now set up the daughter particles

   FillTree(gProduct,gInParam->reaction.product);
   FillTree(gRecoil,gInParam->reaction.recoil); 
}

/////////////////////////////////////////////////////////////////

void Input::FillTree(Particle *theParticle, int itsIndex)
{
// this creates and fills the rest of the particles in the reaction

       	// set pointer to the particle's database info
        theParticle->mInfo = &(gDBPtr->particle[itsIndex]);

        // now do the rest of the initalization for the particle after mInfo
        // is set
        theParticle->InitParticle();

        // tell the various keepers about this particle.

        gMotherKeeper->keep(theParticle);

        if (theParticle->mInfo->returnData) gUserKeeper->keep(theParticle);
	if (theParticle->mInfo->isPhoton) gPhotonKeeper->keep(theParticle);
        
        // set the number of daughters counter

      	theParticle->mNumDaughters = 
              (gDBPtr->particle[itsIndex]).decay[0].nProducts;

        // now make kids

	for (int i=0; i < theParticle->mNumDaughters; i++)
	{

 	// first allocate a new particle in memory

	  Particle *newGuy = NULL;
	  newGuy = new(Particle);
	  assert(newGuy != NULL);

	// now set theParticle's daughter pointer to this newGuy

          theParticle->mDaughter[i] = newGuy;

        // Let the daughter know it's parent

          newGuy->mParent=theParticle;

        // now FillTree for the newGuy

          FillTree(newGuy,gDBPtr->particle[itsIndex].decay[0].product[i]);

	}
}

/////////////////////////////////////////////////////////////////

void Input::WriteHeader(void)
{
}

