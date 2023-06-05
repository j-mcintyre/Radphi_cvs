///////////////////////////////////////////////////////////////////////////////
//                                            
//  Chees.C                                    
//
//  14 June 1995                                 
//                                                
//  This file contains the three main libChees subroutines--
//
//  CutChees : initialize
//  GenerateOne : generates one event
//  ScrubChees : cleans up.  
//   Call before end of program or before reading new input file. 
//   Some parameters can be changed without Scrubbbing (things that don't 
//   depend on histograms).
// 
///////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// header files

#ifndef CHEES_INCLUDED
 #include "Chees.h"
#endif

#include "Input.h"
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include <wabbit.h>
#include <param.h>
}

#include "Particle.h"
#include "Keeper.h"
#include "Kinmaster.h"
#include "util.h"


////////////////////////////////////////////////////////////////////////
// global variables

Keeper *gUserKeeper;
Keeper *gPhotonKeeper;
Keeper *gMotherKeeper;

    /* 
       The Keeper is in charge of holding on to particles. 
       The Keeper also does the smearing, cutting, and
       four-vector output.

       gUserKeeper holds on to (and outputs) user specified particles.
       gPhotonKeeper holds on to (and cuts and smears) the photons.
       gMotherKeeper holds on to (and Lorentz transforms) everyone.

    */

Particle *gBeam;
Particle *gTarget;
Particle *gProduct;
Particle *gRecoil; 

Four      *gLabFrameVec;
    
    /* 
       The beam, target, product and recoil are used in every reaction.
       The rest of the Particles, which depend on the input file,  are 
       created in Initialize().
       
       The labFrame vector is used to transform from the CM frame to the
       lab frame once we're done generating.  It starts out with p=0 and 
       some dummy mass.  (We can't use the gTarget vector to specify
       the lab frame because the target
       particles may have some (Fermi) momentum.
       
    */

Input *gUserInput;

    /* 
       The Input object reads the input file and holds on to the Input
       parameters.
    */

FILE *cprnt;

Kinmaster *gKinmaster;

    /*
       The Kinmaster does the kinematics that are associated with
       groups of particles (as opposed to operations with single
       particles, which are handle in the Particle class).
       E.g., it choose the momentum transfer (t) for each event.
    */


mc_param_t *gInParam;     /* input parameter files */
particles_db_t  *gDBPtr;



/////////////////////////////////////////////////////////////////////
// private function prototype

void CreateGlobalObjects(void);
void Generate(Particle *, Particle *, Particle *, Particle *, 
	      Keeper *, Keeper *);

/////////////////////////////////////////////////////////////////////
// Rob's run number hack - 

static int runNo=0;


void setCheesRunNo(int userRunNo)
{
  runNo = userRunNo;  			  
}		  

/////////////////////////////////////////////////////////////////////
//The routines

/*

CutChees: call once before generating.  It reads in the input file.

You must allocate space for the mc_param_t in your program
 and pass CutChees the pointer.
 
You must also allocate a pointer to a particles_db_t and pass
 a pointer to it to CutChees (i.e. you give it a pointer to a pointer).
 
Your code will look something like this:

 mc_param_t parameters;        
 particles_db_t *particleDBPtr;         

 CutChees(inputFilename, &parameters, &particleDBPtr);
 
CHANGING PARAMTERS DURING RUNTIME: 
Once you read in the input file you can change numerical parameters in memory
during runtime.  E.g., you can change the parent mass.  For other parameters,
such as input histograms, you will have to clean up and reinitialize OR
use a special routine (which you might have to write).  E.g., to change
t-distribution histograms, call ChangeTHisto (in this file).  If you write
such a routine please place it in this file.

*/

void CutChees(char *inFile, mc_param_t *paramPtr, 
                            particles_db_t **dbHandle)
{

  /* first create the input object */
  gUserInput = NULL;
  gUserInput = new(Input);
  assert(gUserInput != NULL);

  gInParam = paramPtr;  /* user had to allocate this storage already;
                           we need to point to this storage.
                         */

  /* now read the input file from inFile */
  if(gUserInput->ReadInput(inFile)) /* non zero if it fails */
    {
      printf("\n\nSorry, I had trouble reading your input file.\n\n");
      printf("Input filename: %s\n\n\n",inFile);
      exit(EXIT_FAILURE);
    }

  *dbHandle = gDBPtr;  /* ReadInput allocated storage and set gDBPtr to it.
                          Now set the user's pointer to point to it as well.
                       */ 

  /* now create the space for all of the other global ojects such as the
     Kinmaster, the reaction particles and the Keepers */
  CreateGlobalObjects();


  /* next we wish to setup the reaction and fill in the decay products */
  gUserInput->SetUpParticles();


  /* now open the file for output (if C or C++ format) */
  if((gInParam->codes.outFileType == OUT_C) ||
     (gInParam->codes.outFileType == OUT_CPLUSPLUS))
    cprnt = fopen(gInParam->codes.outFile, "w");

  if(gInParam->codes.geometry == GEOM_E852)
    {   
      fprintf(stderr,"E852 GEOMETRY is not supported in this release. (Chees.C)\n");
    };
  
  /* finally book the histos if the user has specified */
  if(gInParam->codes.doHisto)
    gPhotonKeeper->InitHistos();

  /* that's all folks! */

}

//////////////////////////////////////////////////////////////////////

/*
   the next routine takes a cheesInfo struct, which the user has created and 
   allocated space for, as an argument, generates one event, and then stores
   a list of 4-vectors for the event in the cheesInfo structure.  The routine
   returns an integer 1 if the event passed the geometry cuts and 0 if it
   doesn't
*/

int GenerateOne(cheesInfo_t *userInfo)
{
  
  int failedCuts = 0;
  int doOutput = 1;


  /* first generate all of the four vectors for the reaction and put them in 
     the keepers.  These 4-vectors are created in the GJ frame */
     
  Generate(gBeam, gTarget, gProduct, gRecoil, 
                                gUserKeeper, gPhotonKeeper);

  /* now transform from the GJ frame to the overall CM frame */
  
  gMotherKeeper->RotateGJtoCM();
    
  /* now transform everyone into the lab frame. */
  /* we don't go to the target CM anymore because Fermi
     momentum may be turned on. */
     
  gMotherKeeper->TransformAll(gLabFrameVec);

  /* now reset the lab frame vector (MotherKeeper doesn't
      know about him) */
      
  gLabFrameVec->Set(0.0,0.0,0.0,1.0);  

  /* do the geometry cuts and plot the accepted events if user wishes */

  /* steal a copy of the unsmeared photon four vectors */
  /* check for possible negative photon energies */    

  userInfo->nGammas= gPhotonKeeper->mNumParts;
  
      
  for(int i = 0; i < userInfo->nGammas; ++i){

   /* first check for neg photon energies; if so, flag it and change E */

   if (gPhotonKeeper->mData[i]->mFour.vec.t < 0.0)
   {
     printf("\nDoh! Found negative photon energy!!\n");
     gPhotonKeeper->mData[i]->mFour.print();

     double magp = gPhotonKeeper->mData[i]->mFour.mag3();
     gPhotonKeeper->mData[i]->mFour.vec.t = magp;

     /* failedCuts |= neg_energy_warn; */ /* WARN THE USER! */

   }; /* end negative energy check */

   userInfo->pRawGammas[i].space.x = gPhotonKeeper->mData[i]->mFour.vec.space.x;
   userInfo->pRawGammas[i].space.y = gPhotonKeeper->mData[i]->mFour.vec.space.y;
   userInfo->pRawGammas[i].space.z = gPhotonKeeper->mData[i]->mFour.vec.space.z;
   userInfo->pRawGammas[i].t = gPhotonKeeper->mData[i]->mFour.vec.t;
  } /* end loop over gammas */

  switch (gInParam->codes.geometry)
  {
  	case GEOM_E852:
  	  	 printf("\n\nSorry, E852 geometry not implemented (Chees.C).\n\n");
      	 exit(EXIT_FAILURE);
  		break;
  		
  	case GEOM_CEBAF_328:
  	case GEOM_CEBAF_620:
	    failedCuts = 0;
	    break;
  		
  	case GEOM_GAMS2000:
  	  	 printf("\n\nSorry, GAMS2000 geometry not implemented (Chees.C).\n\n");
      	 exit(EXIT_FAILURE);
		break;
  		
  	case GEOM_GAMS4000:
  		 printf("\n\nSorry, CEBAF geometry not implemented (Chees.C).\n\n");
      	 exit(EXIT_FAILURE);
      	break;
  		
  	default:
  	  	 printf("\n\nError: I don't know what geometry you've specified.\n");
  	  	 printf("(Chees.C)\n\n");
      	 exit(EXIT_FAILURE);
  		break;
  		
  } /* end geometry switch*/
  
  if(gInParam->codes.doHisto)
    {
      if (!failedCuts) 
	     gPhotonKeeper->FillHistos(ACCEPTED_HISTOS);
      
      gPhotonKeeper->FillHistos(GENERATED_HISTOS);
    }

  
  /* put the end result in userInfo.  As is the case with the writing to a 
     file, if the user has not specified the return all of the 4-vectors, 
     then only the accepted ones will be returned */

      userInfo->nVectors = gUserKeeper->mNumParts;
      userInfo->globalAcceptBits = failedCuts;
      userInfo->principleVertex.x = gProduct->mOrigin.x;
      userInfo->principleVertex.y = gProduct->mOrigin.y;
      userInfo->principleVertex.z = gProduct->mOrigin.z;
      userInfo->missingMassSq = SQR(gRecoil->mMass);
      userInfo->genProductMass = gProduct->mMass;

   /* now for the vectors */

  if(gInParam->codes.returnAll)
    {
      for(int i = 0; i < userInfo->nVectors; ++i)
      {
	userInfo->data[i].p.space.x = gUserKeeper->mData[i]->mFour.vec.space.x;
	userInfo->data[i].p.space.y = gUserKeeper->mData[i]->mFour.vec.space.y;
	userInfo->data[i].p.space.z = gUserKeeper->mData[i]->mFour.vec.space.z;
	userInfo->data[i].p.t = gUserKeeper->mData[i]->mFour.vec.t;

        userInfo->data[i].acceptBits = gUserKeeper->mData[i]->mPersonalAcceptBits;
	userInfo->data[i].id = gUserKeeper->mData[i]->mInfo->id;
	
	/* Now we need to set the "type" field */

	userInfo->data[i].type = 0;
	if(gUserKeeper->mData[i] == gBeam){
	  userInfo->data[i].type |= CHEES_BEAM;
	}
	else if(gUserKeeper->mData[i] == gTarget){
	  userInfo->data[i].type |= CHEES_TARGET;
	}
	else if(gUserKeeper->mData[i] == gRecoil){
	  userInfo->data[i].type |= CHEES_RECOIL;
	}
	else if(gUserKeeper->mData[i] == gProduct){
	  userInfo->data[i].type |= CHEES_PRODUCT;
	}
	
	/* Is this is particle we want in an ESR (i.e. PWA'able) 
	   How this is determined :
	   parent is an "Unknown" i.e. a0,f0,etc
	   OR if product is known, use it's children
	   */
	 

	else if(userInfo->data[i].id != Unknown){

	  /* Parent is an "unknown" - standard case */

	  if(gProduct->mInfo->id == Unknown){
	    if(gUserKeeper->mData[i]->mParent->mInfo->id == Unknown)
	      userInfo->data[i].type |= CHEES_ESRPART;
	  }

	  else{
	    /* The parent is "known" - keep the children */

	    if(gUserKeeper->mData[i]->mParent == gProduct){
	      userInfo->data[i].type |= CHEES_ESRPART;	      
	    }
	  }
	}
	
	if(gUserKeeper->mData[i]->mNumDaughters != 0){
	  userInfo->data[i].type |= CHEES_UNSTABLE;
	}
	else{
	  
	  /* Does not decay, must be in the final state,
	     unless it was beam or target */

	  if((gUserKeeper->mData[i] != gBeam) &&
	     (gUserKeeper->mData[i] != gTarget))
	    userInfo->data[i].type |= CHEES_FSP;
	}
	
      }
    }
  else if(!failedCuts)
    {
      for(int i = 0; i < userInfo->nVectors; ++i)
      {	
	userInfo->data[i].p.space.x = gUserKeeper->mData[i]->mFour.vec.space.x;
	userInfo->data[i].p.space.y = gUserKeeper->mData[i]->mFour.vec.space.y;
	userInfo->data[i].p.space.z = gUserKeeper->mData[i]->mFour.vec.space.z;
	userInfo->data[i].p.t = gUserKeeper->mData[i]->mFour.vec.t;
	userInfo->data[i].id = gUserKeeper->mData[i]->mInfo->id;
      };
    }
  

  /* now write the four-vectors out to a file if the user 
     wants that kind of action */
  switch(gInParam->codes.outFileType)
    {
    case OUT_NONE:
      break;
    case OUT_C:
      if(gInParam->codes.returnAll)
	gUserKeeper->cprint(cprnt, failedCuts);
      else if(!failedCuts)
	gUserKeeper->cprint(cprnt, failedCuts);
      break;
    case OUT_CPLUSPLUS:
      if(gInParam->codes.returnAll)
	gUserKeeper->cprint(cprnt, failedCuts);
      else if(!failedCuts)
	gUserKeeper->cprint(cprnt, failedCuts);
      break;
    case OUT_FORTRAN:
      if(doOutput)
	{
	  fprintf(stderr, "Sorry, the output format you are requesting ");
	  fprintf(stderr, "has not currently been implemented.\n");
	  fprintf(stderr, "No output will be produced.\n\n");
	  doOutput = 0;
	}
      break;
    case OUT_ASCII:
      if(doOutput)
	{
	  fprintf(stderr, "Sorry, the output format you are requesting ");
	  fprintf(stderr, "has not currently been implemented.\n");
	  fprintf(stderr, "No output will be produced.\n\n");
	  doOutput = 0;
	}
      break;
    default:
      fprintf(stderr, "Incorrect output format specification\n");
      fprintf(stderr, "Check your input file for errors\n");
      fprintf(stderr, "Exiting...\n\n");
      exit(EXIT_FAILURE);
    }
  
  
  /* finally return the global cuts bitmask (the logical OR of
     all the individual cutmasks.  Zero means everything is OK. */
  return (failedCuts);

}
  
////////////////////////////////////////////////////////////////////////

/*
   The ScrubChees routine will get rid of all of the memory it originally 
   allocated for the global objects.  If the user changes the parameters and 
   wishes to reinitialize, then he/she must call ScrubChees to deallocate the 
   memory before initializing.  Otherwise, memory leaks will occur.
   ScrubChees will also write diagnostic histos to a bbook file if the user 
   has requested so, and it will run efud on that file to make an hbook file 
   out of it.  The output is put in the file histos.bbook and histos.hbook
*/

void ScrubChees(void)
{
  
  /* nuke the global objects that have been created */

  // Keepers

  delete(gUserKeeper);
  delete(gPhotonKeeper);
  delete(gMotherKeeper);

  // Particles

  delete(gBeam);
  delete(gTarget);
  delete(gProduct);
  delete(gRecoil);

  // Input

  delete(gUserInput);

  // Kinmaster

  delete(gKinmaster);

  // Close output files (if necessary).

  if((gInParam->codes.outFileType == OUT_C) ||
     (gInParam->codes.outFileType == OUT_CPLUSPLUS))
    fclose(cprnt);

  /* now write the histos to a file and run efud if desired */
  if(gInParam->codes.doHisto)
    {
      gPhotonKeeper->WriteHistos();
      runEfud("histos.bbook");
    }

}

//////////////////////////////////////////////////////////////////////

/* 
   This puppy lets the user change t-distribution histograms while running,
   without having to ScrubChees.
*/


void ChangeTHisto(int histoID,char *filename)
{
	if (gKinmaster->mUserTPicker->mData) 
          free(gKinmaster->mUserTPicker->mData);  /* nuke old one */
    
        gKinmaster->mUserTPicker->InitNosPicker(histoID,filename);
}


////////////////////////////////////////////////////////////////////////
/*
   This function allocates space for the global C++ objects
*/

void CreateGlobalObjects(void)
{

// create the keepers

  gUserKeeper = NULL;
  gUserKeeper = new(Keeper);
  assert(gUserKeeper != NULL);
        
  gPhotonKeeper = NULL;
  gPhotonKeeper = new(Keeper);
  assert(gPhotonKeeper != NULL);
 
  gMotherKeeper = NULL;
  gMotherKeeper = new(Keeper);
  assert(gMotherKeeper != NULL);

// create the standard particles

  gBeam = NULL;
  gBeam = new(Particle);
  assert(gBeam != NULL);

  gTarget = NULL;
  gTarget = new(Particle);
  assert(gTarget != NULL);
  
  gProduct = NULL;
  gProduct = new(Particle);
  assert(gProduct != NULL);

  gRecoil = NULL;
  gRecoil = new(Particle);
  assert(gRecoil != NULL);
  
  gLabFrameVec = NULL;
  gLabFrameVec = new(Four);
  assert(gLabFrameVec != NULL);
  gLabFrameVec->Set(0.0,0.0,0.0,1.0 /*p=0, dummy mass */);
  
// create the Kinmaster 

  gKinmaster = NULL;
  gKinmaster = new(Kinmaster);
  assert(gKinmaster != NULL);
}





