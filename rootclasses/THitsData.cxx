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
//   from TTree THitsData/Radphi
//   found on file: simData1.root and subsequently modified.
//   It provides an interface to the hit data section of 
//   a data tree.  One way to get such a tree is to start
//   with an itape file, use mctuple or ntuple as appropriate
//   and get a paw ntuple.  Then use h2root to convert the 
//   paw ntuple to a ROOT TTree.
//   This class contains room enough for any 1 event, and
//   through the GetEvent() member, it can load any event
//   from a TTree.
//   It contains any histogram making, data displaying, or
//   calculation functions that use the hit data.
//////////////////////////////////////////////////////////

#include <iostream.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include "THitsData.h"

ClassImp(THitsData)

THitsData::THitsData(TTree *tree)
{
// "marries" this THitsData to tree, and sets the branches of tree to 
// point to the data members of this
// by default, it will search for the tree h10 in the file simData.root
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("simData.root");
      if (!f) {
         f = new TFile("simData.root");
      }
      tree = (TTree*)gDirectory->Get("h10");

   }
   Init(tree);
}

Int_t THitsData::GetEvent(Int_t event)
{
// Read specified event from the Tree into data members
   if (fTree) return fTree->GetEvent(event);
   else       return 0;
}

void THitsData::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fTree = tree;
      fTree->SetBranchAddress("runno",&runno);
      fTree->SetBranchAddress("eventno",&eventno);
      fTree->SetBranchAddress("ismc",&ismc);
      fTree->SetBranchAddress("nhbsd",&nhbsd);
      fTree->SetBranchAddress("chbsd",chbsd);
      fTree->SetBranchAddress("Ebsd",Ebsd);
      fTree->SetBranchAddress("ntbsd",ntbsd);
      fTree->SetBranchAddress("t1bsd",t1bsd);
      fTree->SetBranchAddress("nhbgv",&nhbgv);
      fTree->SetBranchAddress("chbgv",chbgv);
      fTree->SetBranchAddress("Ebgvdwn",Ebgvdwn);
      fTree->SetBranchAddress("Ebgvup",Ebgvup);
      fTree->SetBranchAddress("tbgvdwn",tbgvdwn);
      fTree->SetBranchAddress("tbgvup",tbgvup);
      fTree->SetBranchAddress("nhcpv",&nhcpv);
      fTree->SetBranchAddress("chcpv",chcpv);
      fTree->SetBranchAddress("Ecpv",Ecpv);
      fTree->SetBranchAddress("ntcpv",ntcpv);
      fTree->SetBranchAddress("t1cpv",t1cpv);
      fTree->SetBranchAddress("nhupv",&nhupv);
      fTree->SetBranchAddress("chupv",chupv);
      fTree->SetBranchAddress("Eupv",Eupv);
      fTree->SetBranchAddress("ntupv",ntupv);
      fTree->SetBranchAddress("t1upv",t1upv);
      fTree->SetBranchAddress("nhlgd",&nhlgd);
      fTree->SetBranchAddress("chlgd",chlgd);
      fTree->SetBranchAddress("Elgd",Elgd);
      fTree->SetBranchAddress("ntimes",&ntimes);
      fTree->SetBranchAddress("le",le);
}

void THitsData::Show(Int_t event) const
{
// Print contents of event.
// If event is not specified, print current event
   if (!fTree) return;
   fTree->Show(event);
}

TH1 *THitsData::Plotnhbgv()
{
// returns a TH1 of the number of upstream hits in the BGV over threshold in 
// each event. Typical usage in ROOT interactive would be hits.Plotnhbgv()->Draw()
  cout << "Plotting number of upstream hits in the BGV over threshold per event." << endl;
  TMath math;
  TH1D *hist=new TH1D ("nhbgv","Upstream hits in the BGV over threshold per event.",10,0,10);
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("nhbgv",1);
  fTree->SetBranchStatus("tbgvdwn",1);
  fTree->SetBranchStatus("tbgvup",1);
  fTree->SetBranchStatus("chbgv",1);
  Int_t numEvents=fTree->GetEntries();
  cout << numEvents << " events to be processed." << endl;
  Int_t hitNum;
  Int_t numGoodHits;
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    GetEvent(eventNum);
    numGoodHits=0;
    for(hitNum=0;hitNum<nhbgv;hitNum++) {
      if (tbgvdwn[hitNum]<1e-6&&tbgvup[hitNum]<1e-6) numGoodHits++;
    }
    hist->Fill(numGoodHits);
    if (eventNum%1000==0)
      cout << eventNum << " events processed.\r" << flush;
  }
  cout << endl << "Done." << endl;
  return hist;
}
