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
//   from TTree TReconstruction/Radphi
//   found on file: simData1.root and modified signifigantly
//   This class provides an interface to the reconstruction 
//   portion of an event data 
//   tree, which should contain data for number of events.
//   Such a tree can be constructed from a itape file by
//   first using ntuple or mctuple to generate a paw
//   ntuple, and then using h2root to generate a ROOT TTree
//   This class contains room for 1 event from this tree,
//   but can load any event from the TTree using the GetEvent()
//   member.
//   It contains any functions for displaying or calculating
//   from the reconstrucion data.
//   note: uses a +--- metric.
//////////////////////////////////////////////////////////

#include <iostream.h>
#include <ieeefp.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TSortedList.h>

#include "TChannel.h"
#include "TReconstruction.h"
#include "TResolver.h"
#include "TTable.h"

#define PI 3.1415926535897932

ClassImp(TReconstruction)

void TReconstruction::FindPi0(Int_t eventNum)
{
  // searches all pairs of photons in pvect of event eventNum,
  // looking for pi0s, and prints the results
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("pvect",1);
  fTree->SetBranchStatus("nphot",1);
  GetEvent(eventNum);
  for(Int_t photon1=0;photon1<nphot;photon1++) {
    TChannel gamma1(kGamma1,pvect[photon1]);
    for(Int_t photon2=photon1+1;photon2<nphot;photon2++) {
      TChannel gamma2(kGamma1,pvect[photon2]);
      TChannel pi0(kPi02,gamma1,gamma2);
      if (pi0.IsInRange()) {
        cout << "got a pion on event " << eventNum << endl;
	cout << "the mass found is " << pi0.M() << endl;
      }
    }
  }
}

void TReconstruction::FindPi0()
{
  // runs FindPi0(EChannelKind) for every event.  
  // note: Prints a lot of stuff
  Int_t numEntries;
  numEntries=fTree->GetEntries();
  cout << numEntries << " to be searched for pi0s" << endl;
  for(Int_t entryNum=0;entryNum<numEntries;entryNum++) {
    FindPi0(entryNum);
  }
}

TH1 *TReconstruction::FitParticle(EChannelKind particle)
{
// will fit a gaussian around the mass of particle to a plot
// of the masses of pairs of photons.
// note that if you set particle to anything that doesn't
// decay into 2 photons, you'll get nonsense
  Double_t mass=TChannel::M0(particle);
  Double_t sigma=TChannel::Sigma(particle);
  TH1 *hist;
  hist = (TH1D*)gDirectory->Get("Mof2");
  if (hist==NULL) 
    hist=PlotMof2();
  TF1 f1("massFit","gaus(0)+[3]",TChannel::LowerBound(particle),
				 TChannel::UpperBound(particle));
  f1.SetParNames("amplitude","mean","sigma","pedestal");
  f1.SetParameters(500,mass,sigma,50);
  hist->Fit("massFit","NR");
  return(hist);
}

TH1 **TReconstruction::FindSigma(Int_t maxEvents)
{
  // creates a mass plot for each channel
  // and fits each to a constant + gaussian,
  // and creates a table of results
  // reads at most maxEvents from its tree
  // returns an array of pointers to these histograms
  TH1D **hist=new TH1D* [kNumChannels];
  {
    for(Int_t histNum=0;histNum<kNumChannels;histNum++) {
      TString name,title;
      name=TChannel::Name((EChannelKind)histNum);
      title="The mass spectrum of ";
      title+=name;
      name+="hist";
      hist[histNum]= new TH1D(name.Data(),title.Data(),100,
	TChannel::LowerBound((EChannelKind)histNum),
        TChannel::UpperBound((EChannelKind)histNum));
    }
  }
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("nphot",1);
  fTree->SetBranchStatus("pvect",1);
  Int_t numEvents=fTree->GetEntries();
  numEvents = (numEvents < maxEvents) ? numEvents : maxEvents;
  cout << numEvents << " events to be processed." << endl;
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    GetEvent(eventNum);
    TResolver resolver(this);
    resolver.FindAll();
    for(Int_t channelNum=0;channelNum<kNumChannels;channelNum++) {
      TIter iter(resolver.GetChannelList((EChannelKind)channelNum));
      TChannel *channel;
      while(channel= (TChannel*)iter()) {
        hist[channelNum]->Fill(channel->M());
      }
    }
    if (eventNum%100==0)
      cout << eventNum << " events processed.\r" << flush;
  }
  cout << endl;
  Double_t fitResults[kNumChannels][4];
  TF1 f("gwp","[0]+gaus(1)");
  f.SetParNames("pedistal","amplitude","mass","sigma");
  cout << "Fitting data...            " << endl;
  for(Int_t channelNum=0;channelNum<kNumChannels;channelNum++) {
    f.SetRange(TChannel::LowerBound((EChannelKind)channelNum),
	       TChannel::UpperBound((EChannelKind)channelNum));
    f.SetParameters(50,100,TChannel::M0((EChannelKind)channelNum),
  			   TChannel::Sigma((EChannelKind)channelNum));
    hist[channelNum]->Fit("gwp","RN");
    f.GetParameters(fitResults[channelNum]);
  }
  Int_t cw[]={10,8,8,8,8,8,6};
  TTable table(7,cw);
  table.PrintBar();
  table.Print("Name");
  table.Print("Pedestal");
  table.Print("Amp.");
  table.Print("M-M0");
  table.Print("Sigma");
  table.Print("sig/bg");
  table.Print("number");
  table.PrintBar();
  TMath math;
  for(channelNum=0;channelNum<kNumChannels;channelNum++) {
    if (hist[channelNum]->GetEntries()==0) continue;
    table.Print(TChannel::Name((EChannelKind)channelNum));
    table.Print(fitResults[channelNum][0]);
    table.Print(fitResults[channelNum][1]);
    table.Print(fitResults[channelNum][2]-
	        TChannel::M0((EChannelKind)channelNum));
    table.Print(fitResults[channelNum][3]);
    table.Print(fitResults[channelNum][1]/fitResults[channelNum][0]);
    table.Print((Int_t)hist[channelNum]->GetEntries());
  }
  table.PrintBar();
  return((TH1 **)hist);
}

TReconstruction::TReconstruction(TTree *tree)
{
  // "marries" this to tree, and sets tree's reconstruction
  // branches to point to this
  // by default, this look for the tree h10 in the
  // file simData.root
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("simData3.root");
      if (!f) {
         f = new TFile("simData3.root");
      }
      tree = (TTree*)gDirectory->Get("h10");

   }
   Init(tree);
}

Int_t TReconstruction::GetEvent(Int_t event)
{
// Read specified event from the Tree into data members
   if (fTree) return fTree->GetEvent(event);
   else       return 0;
}

void TReconstruction::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fTree = tree;
      fTree->SetBranchAddress("nphot",&nphot);
      fTree->SetBranchAddress("nfrwd",&nfrwd);
      fTree->SetBranchAddress("pvect",pvect);
}

void TReconstruction::Show(Int_t event) const
{
// Print contents of event.
// If event is not specified, print current event
   if (!fTree) return;
   fTree->Show(event);
}

TH1 *TReconstruction::PlotMof2()
{
  // creates a mass plot of every pair of 2 photons.
  cout << "Plotting the mass of every pair of photons." << endl;
  TMath math;
  TH1D *hist=new TH1D("Mof2","The mass of every possible photon pair",100,0,0.8);
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("nphot",1);
  fTree->SetBranchStatus("pvect",1);
  Int_t numEvents=fTree->GetEntries();
  Int_t photonNum1,photonNum2;
  cout << "creating a plot of possible decay parent of 2 photonNum pairs." << endl;
  cout << numEvents << " events total." << endl;
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    GetEvent(eventNum);
    for(photonNum1=0;photonNum1<nphot;photonNum1++) {
      TChannel gamma1(kGamma1,pvect[photonNum1]);
      for(photonNum2=photonNum1+1;photonNum2<nphot;photonNum2++) {
        TChannel gamma2(kGamma1,pvect[photonNum2]);
        TChannel pi0(kPi02,gamma1,gamma2);
        hist->Fill(pi0.M());
      }
    }
  }
  return hist;
}

TH2 *TReconstruction::Plot2Particles()
{
  // creates a 2-d histogram of all possible pairs of pairs of photons
  cout << "Building a 2-d histogram of possible 2 particle states decaying from 2 photons." 
       << endl;
  Int_t p1a,p1b,p2a,p2b; // p1a and p1b make a pair, as do p2a and p2b
  		// p1a<p1b to avoid double counting, likewise p2a<p2b
		// p2a and p2b don't equal either p1a or p1b
  TH2D *hist=new TH2D("2Particles",
 		      "A plot of the masses of every pair of photon pairs",
		      100,0,0.8,100,0,0.8);
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("nphot",1);
  fTree->SetBranchStatus("pvect",1);
  Int_t numEvents=fTree->GetEntries();
  numEvents=(numEvents>10000) ? 10000 : numEvents;
  cout << numEvents << " events to be processed." << endl;
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    GetEvent(eventNum);
    if (!eventNum%100) cout << eventNum << '\r';
    if (nphot>=4) {
      for(p1a=0;p1a<nphot-1;p1a++) {
        TChannel gamma1a(kGamma1,pvect[p1a]);
        for(p1b=p1a+1;p1b<nphot;p1b++) {
          TChannel gamma1b(kGamma1,pvect[p1b]);
	  TChannel pi01(kPi02,gamma1a,gamma1b);
          for(p2a=0;p2a<nphot-1;p2a++) {
	    if ((p2a!=p1a)&&(p2a!=p1b)) {
              TChannel gamma2a(kGamma1,pvect[p2a]);
              for(p2b=p2a+1;p2b<nphot;p2b++) {
		if ((p2b!=p1a)&&(p2b!=p1b)) {
                  TChannel gamma2b(kGamma1,pvect[p2b]);
		  TChannel pi02(kPi02,gamma2a,gamma2b);
	          hist->Fill(pi01.M(),pi02.M(),1.0);
		} //p2b conditional
 	      } //p2b loop
	    } //p2a conditional
	  } //p2a loop
        } //p1b loop
      } //p1a loop
    } //nphot>=4 if statement
  } //event loop
  return(hist);
}

TH2 *TReconstruction::PlotThetaPhiOfPi0()
{
  // makes a 2-d histogram of the theta and phi for
  // the photons of pi0s in the pi0's reference frame
  cout << "Plotting the angular distribution of pi0 decays" << endl;
  TH2D *hist=new TH2D("TPpi0","Theta & phi of decay photons in the CM frame of pi0",
		      100,0,PI,100,0,2*PI);
  Int_t particleNum1,particleNum2;
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("nphot",1);
  fTree->SetBranchStatus("pvect",1);
  Int_t numEvents=fTree->GetEntries();
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    GetEvent(eventNum);
    for(particleNum1=0;particleNum1<nphot-1;particleNum1++) {
      TChannel gamma1(kGamma1,pvect[particleNum1]);
      for(particleNum2=particleNum1+1;particleNum2<nphot;particleNum2++) {
        TChannel gamma2(kGamma1,pvect[particleNum2]);
        TChannel pi0(kPi02,gamma1,gamma2);
        if (pi0.IsInRange()) {
          gamma1.Boost(pi0.BoostVector());
          hist->Fill(gamma1.Theta(),gamma1.Phi(),1.0);
        }
      }
    }
    if (eventNum%1000==0)
      cout << eventNum << " events processed.\r" << flush;
  }
  cout << endl << "done." << endl;
  return hist;
}

TH1 *TReconstruction::PlotChannelMass(EChannelKind kind)
{
  // creates a mass plot of a given channel
  // note: is about as slow as FindSigma(), which does all channels
  TH1D *hist = new TH1D("ChMhist","Mass of the parent of a given channel",
                        100,TChannel::LowerBound(kind),TChannel::UpperBound(kind));
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("nphot",1);
  fTree->SetBranchStatus("pvect",1);
  Int_t numEvents=fTree->GetEntries();
  cout << numEvents << " events to be plotted." << endl;
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    GetEvent(eventNum);
    TResolver resolve(this);
    resolve.FindAll();
    TIter iter(resolve.GetChannelList(kind));
    TChannel *channel;
    while (channel= (TChannel*)iter()) {
      hist->Fill(channel->M());
    }
    if (eventNum%1000==0)
      cout << eventNum << " events processed." << '\r' << flush;
  }
  return(hist);
}

TH1 **TReconstruction::PlotEBGV_ELGD()
{
  // returns an array of 2 pointers to histograms for the LGD and BGV 
  // LGD is hist[0] and BGV is hist[1]
  fTree->SetBranchStatus("*",0);  // disable all branches
  fTree->SetBranchStatus("nphot",1);  // activate branchname
  fTree->SetBranchStatus("nfrwd",1);  // activate branchname
  fTree->SetBranchStatus("pvect",1);  // activate branchname
  TH1D **hist= new TH1D*[2];
  hist[0]= new TH1D("ELGD","forward cluster energies",100,0.01,5);
  hist[1]= new TH1D("EBGV","BGV cluster energies",100,0.01,1);
  if (fTree == 0) return NULL;
  Int_t nentries = fTree->GetEntries();
  Int_t nbytes = 0, nb = 0;
  Int_t clusterNum;
  for (Int_t i=0; i<nentries;i++) {
    nb = fTree->GetEvent(i);   nbytes += nb;
    for(clusterNum=0;clusterNum<nfrwd;clusterNum++) {
      hist[0]->Fill(pvect[clusterNum][0]);
    }
    for(clusterNum=nfrwd;clusterNum<nphot;clusterNum++) {
      hist[1]->Fill(pvect[clusterNum][0]);
    }
  }
  return (TH1**)hist;
}

TH1 *TReconstruction::PlotMofBGVphoton()
{
  // creates a histogram of the mass of the photons in the BGV
  TH1C *hist=new TH1C("MofBGVphotons",
       "Mass of the photons in the BGV, seriously",50,0.0,0.3);
  hist->SetMaximum(10);
  fTree->SetBranchStatus("*");
  fTree->SetBranchStatus("nphot",1);
  fTree->SetBranchStatus("pvect",1);
  fTree->SetBranchStatus("nfrwd",1);
  Int_t numEvents=fTree->GetEntries();
  cout << numEvents << " to be processed" << endl;
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    GetEvent(eventNum);
    for(Int_t photNum=nfrwd;photNum<nphot;photNum++) {
      TChannel gamma(kGamma1,pvect[photNum]);
      hist->Fill(gamma.M());
    }
    if (eventNum%1000==0)
      cout << eventNum << " events processed\r" << flush;
  }
  return hist;
}

TH1 *TReconstruction::PlotThetaOfGammas()
{
  // plots the theta angle of the photons
  TH1D *hist=new TH1D("gammaTheta","Theta of photons",100,0,2);
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("pvect",1);
  fTree->SetBranchStatus("nphot",1);
  fTree->SetBranchStatus("nfrwd",1);
  Int_t numEvents=fTree->GetEntries();
  cout << numEvents << " to be processed" << endl;
  for(Int_t eventNum=0;eventNum<numEvents;eventNum++) {
    if (eventNum%1000==0)
      cout << eventNum << " events processed.\r" << flush;
    fTree->GetEvent(eventNum);
    for(Int_t photNum=0;photNum<nphot;photNum++) {
      Double_t theta=TMath::ACos(pvect[photNum][3]/pvect[photNum][0]);
      if (isnand(theta)) {
        printf("theta=NAN, pz=%e, E=%e\n",pvect[photNum][3],pvect[photNum][0]);
      }
      else
        hist->Fill(theta);
    }
  }
  cout << endl;
  return hist;
}
