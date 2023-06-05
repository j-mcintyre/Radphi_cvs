// collide.C
//
// subroutine used in Bnlgams which generates beam + target --> product +
//                                                                 recoil 
//
// PCF 19 June 1995
//////////////////////////////////////////////////////////////////

#include <math.h>
#include "Particle.h"
#include "Keeper.h"
#include "Kinmaster.h"
#include "util.h"

const double Pi = 3.14159;

// links to globals //////////////////////////////////////////////

extern Kinmaster *gKinmaster;
extern Four      *gLabFrameVec;
extern Keeper    *gMotherKeeper;

// prototypes ////////////////////////////////////////////////////

void collide(Particle*,Particle*,Particle*,Particle*);
double unfRand(void);

//////////////////////////////////////////////////////////////////

// collide //

/*
   This routine generates beam + target --> product + recoil in
   the overall CM frame.  Assumes particles are in lab frame to begin
   with and beam + target are already specified.  Fills in product
   and recoil particles.

*/

void collide( Particle *beam
             ,Particle *target
             ,Particle *product
             ,Particle *recoil
            )

{

int numberOfTries;
  numberOfTries = 0;  // used in check to see if there is enough CM energy
static Four theCmVec;
double w; // total CM energy

do{     
      // Choose masses according to user input specifications.

        beam->ChooseMass();
        target->ChooseMass();
        product->ChooseMass();
        recoil->ChooseMass();

      // Now have the Kinmaster prepare the beam
      // and target momenta. 
       
		gKinmaster->PrepareBeam(beam);
		gKinmaster->PrepareTarget(target);
	    
	  // At this point, the beam and target four vectors
 	  // are completely specified in the lab frame.
 	  
 	  // The "labFrame" 4 vec  starts out as p=0 (rest), E=1.0 (dummy value) --
 	  //  it specifies the lab frame.
 	  // Now take it into the CM frame.
 	  // We need to do this because after all the 4 vecs are generated,
 	  //  they'll be in the CM frame.  The Lorentz x-form determine s
 	  //  the beta for the CM-->lab x-form
 	  //  by looking at labFrame vector (gLabFrameVec). 
 	  //  The actual transform happens in Chees.C. 
 	   
 	   	
 	   theCmVec = beam->mFour + target->mFour;
  	   
  	   // compute and store total CM energy
  	  
  	  	 w = theCmVec.mag4();
	     gKinmaster->mTotCMEnergy = w;
	     
	  numberOfTries++;
	  if (numberOfTries % 1000 == 0)
	  { 
	    printf("\n I'm having trouble getting enough CM energy to make your product.\n");
        printf("Check your input file... beam momentum high enough?\n\n");
        numberOfTries = 0;
      }
 }/* end do loop*/	     
	   // make sure we have enough CM energy
	       while (w < (product->mMass + recoil->mMass) ) ;
	    
 	    gLabFrameVec->Set(0.0,0.0,0.0,1.0 /* dummy mass*/);
		gKinmaster->LorentzXForm1( &theCmVec,gLabFrameVec );
		
	  // take beam and target into the CM frame.
	  
		gKinmaster->LorentzXForm1( &theCmVec,&(beam->mFour) );
		gKinmaster->LorentzXForm1( &theCmVec,&(target->mFour) );

 	  // Set up GJ --> CM rotation matrix.
 	  //
      // We rotate GJ --> CM 
 	  // after everything is generated (this happens in Chees.C)
  	  // and sitting in the GJ frame.
  	  // 
 	  // The rotation is first about the y axis (eliminates x
 	  //  component of the beam) by angle thetaXZ,
 	  //  then about the x axis (eliminames y component) by
 	  //  angle thetaYZ.  (The inverse rotation we're constructing now
 	  //  replaces the components in reverse order.)
  	    
 	
 	    
 	    double bx,by,bz;
 	     bx = beam->mFour.vec.space.x;
 	     by = beam->mFour.vec.space.y;  // short-hand for formulas below.
 	     bz = beam->mFour.vec.space.z;
 	    
	    double costhetaXZsq = SQR(bz) / (SQR(bx) + SQR(bz)); // actually cos^2        
 	    double sinthetaXZsq = 1.0 - costhetaXZsq;
 	       double   sinthetaXZ = sqrt(sinthetaXZsq);
 	       double   costhetaXZ = sqrt(costhetaXZsq);
 	          
 	     
 	    double zprime = bx * sinthetaXZ + bz * costhetaXZ;

 	    			   // zprime is new beam z momentum after rotation
 	    			   // about the y axis.
 	    			   
   	    double costhetaYZsq = SQR(zprime) / (SQR(by) + SQR(zprime)); // actually cos^2T
 		double sinthetaYZsq = 1.0 - costhetaYZsq;
              double costhetaYZ = sqrt(costhetaYZsq);
 	          double sinthetaYZ = sqrt(sinthetaYZsq);
 
 	          
 	    // Get the signs right.
 	    // Positive angle rotations go from the x or y axis towards the
 	    // z axis.  The sign makes a difference in the sine terms.
 	    
 	          
 	    sinthetaXZ = ( beam->mFour.vec.space.x > 0.0 ? 
 	                   ABS(sinthetaXZ) :
 	                  -ABS(sinthetaXZ) );

 		sinthetaYZ = ( beam->mFour.vec.space.y > 0.0 ? 
 	                   ABS(sinthetaYZ) :
 	                  -ABS(sinthetaYZ) );
 	        
 	    	   
 	  // Now that we have the rotation angles, tell the MotherKeeper
 	  // to set up the rotation matrix that takes GJ vecs to overall CM vecs.
 	  
 	  gMotherKeeper->
 	     ConstructRotMatrix(costhetaXZ,costhetaYZ,sinthetaXZ,sinthetaYZ);
 	  
 	  // Now set up the matrix that takes our CM vectors to the GJ frame.
 	  // We only have to rotate beam and target (other vecs haven't
 	  // been created yet).
 	  
 	  static SqMatrix3 CMtoGJRot; // aligns CM vecs with the z-axis (GJ frame)
 	  
 	 Three row1( costhetaXZ,                        0, 	           -sinthetaXZ      );
 	 Three row2( -sinthetaYZ*sinthetaXZ,	   costhetaYZ,	-sinthetaYZ*costhetaXZ  );
 	 Three row3(  sinthetaXZ*costhetaYZ,       sinthetaYZ,	 costhetaXZ*costhetaYZ  );
 	  
 	  CMtoGJRot.Set(row1,row2,row3); // fills in matrix
  	  
 	  beam->mFour = CMtoGJRot * (beam->mFour);     // matrix * vector (cool eh?)
 	  target->mFour = CMtoGJRot * (target->mFour);
 	   

	  // The following used to set up beam, target in the GJ frame.
	  // If the Lorentz x-forms had no systematic error, we
	  // could do it this way to save time.  Instead we've
	  // boosted and rotated by hand.
      //  gKinmaster->
      //    FillCMVectors(w,1.0,0.0,beam,target);
  
 
 
      // Set the reaction vertex (filled in by PrepareBeam 
      //  in the beam's mOrigin vector).
      //  Propagation of the vertex for stable particles
      //  is handled in the various Particle::decayx() functions.

        product->mOrigin.x = beam->mOrigin.x;
        product->mOrigin.y = beam->mOrigin.y;
        product->mOrigin.z = beam->mOrigin.z;

        recoil->mOrigin.x = beam->mOrigin.x;
        recoil->mOrigin.y = beam->mOrigin.y;
        recoil->mOrigin.z = beam->mOrigin.z;
 
 
      // Now we need to fill the product and recoil vectors.
      // First find t and use it to get cosTheta between beam
      // and product.  Pick phi uniformly.

double prodCMEnergy = gKinmaster->GetProductCMEnergy(w,product,recoil);
double prodCMMomentum = sqrt(SQR(prodCMEnergy) - SQR(product->mMass));   

double t = gKinmaster->
            ChooseT(prodCMEnergy,prodCMMomentum,beam,product);

double  cosTheta = 
 (t - SQR(beam->mMass) - SQR(product->mMass) + 2*(beam->Energy())*prodCMEnergy)/
                  (2*(beam->Momentum())*prodCMMomentum);

double phi = 2.0 * Pi * unfRand();
        
      // Use this info to fill the product and recoil four vectors.

        gKinmaster->
          FillCMVectors(w,cosTheta,phi,product,recoil);

}

//////////////////////////////////////////////////////////////////


        
                    






