//  DAW

#include "Keeper.h"
#include "Particle.h"


/// prototypes ///

double collide(Particle*,Particle*,Particle*,Particle*);

////////////////////////////////////////////////////


void Generate(Particle *beam
             ,Particle *target
             ,Particle *product
             ,Particle *recoil

             ,Keeper *userKeeper
             ,Keeper *photonKeeper)
{

   // collide everyone who's collidin'

  collide(beam, target, product, recoil);

   // decay everyone who's decayin'

  if ( product->mInfo->nDecays )
    product->Decay(userKeeper, photonKeeper,1);  /* the 1 means
                                                    do GJ Inverse stuff */

  if ( recoil->mInfo->nDecays )  
    recoil->Decay(userKeeper, photonKeeper,0);

}

