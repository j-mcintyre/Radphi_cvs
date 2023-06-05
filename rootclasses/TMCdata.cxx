//---Author: Tom Bogue
//---Version: 1.00/00

/*************************************************************************
 * Copyright(c) 1999, University of Connecticut, All rights reserved.    *
 * Author: Tom Bogue                                                     *
 *                                                                       *
 * Permission to use, copy, modify and distribute this software and its  *
 * documentation for non-commercial purposes is hereby granted without   *
 * fee, provided that the above copyright notice appears in all copies   *
 * and that both the copyright notice and this permission notice appear  *
 * in the supporting documentation. The author makes no claims about the *
 * suitability of this software for any purpose.                         *
 * It is provided "as is" without express or implied warranty.           *
 *************************************************************************/

//////////////////////////////////////////////////////////
//   This class provides an interface to the MC data 
//   portion of a tree created by h2root from a paw nutple
//   that is created by mctuple.  To go over that again,
//   start w/ a itape file from GEANT.  use mctuple
//   to convert this to a paw ntuple.  Then use h2root
//   to convert this to a root TTree.  This class was 
//   origionally produced by using the MakeClass() function
//   member of that tree.  It was modified into this.  It
//   should contain any ploting functions for the MC data
//   as well as any functions to calculate or display the
//   MC data.
//   Note that this class has the capacity to store only
//   1 event, but w/ the GetEvent() function, you may
//   load any given event from the TTree.
//////////////////////////////////////////////////////////

#include <iostream.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include "TMCdata.h"

ClassImp(TMCdata)

TMCdata::TMCdata(TTree *tree)
{
// constructs the TMCdata.  Since a TMCdata is linked for "life" to a TTree
// you should give this a TTree *.  If you use the default value for tree,
// it will try to find the TTree on its own, using a default filename of simData.root
// and a default TTree name of h10
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("simData.root");
      if (!f) {
         f = new TFile("simData.root");
      }
      tree = (TTree*)gDirectory->Get("h10");

   }
   Init(tree);
}

Int_t TMCdata::GetEvent(Int_t event)
{
// Read specified event from the Tree into data members
   if (fTree) return fTree->GetEvent(event);
   else       return 0;
}

void TMCdata::Init(TTree *tree)
{
//   Set branch addresses
// i.e., sets the given branches on fTree to point to this class
   if (tree == 0) return;
   fTree = tree;
      fTree->SetBranchAddress("ni",&ni);
      fTree->SetBranchAddress("type",&type);
      fTree->SetBranchAddress("wgt",&wgt);
      fTree->SetBranchAddress("kindi",kindi);
      fTree->SetBranchAddress("momi",momi);
      fTree->SetBranchAddress("nv",&nv);
      fTree->SetBranchAddress("kindv",kindv);
      fTree->SetBranchAddress("origv",origv);
      fTree->SetBranchAddress("nf",&nf);
      fTree->SetBranchAddress("kindf",kindf);
      fTree->SetBranchAddress("vertf",vertf);
      fTree->SetBranchAddress("momf",momf);
      
}

void TMCdata::Show(Int_t event) const
{
// Print contents of event.
// If event is not specified, print current event
   if (!fTree) return;
   fTree->Show(event);
}
