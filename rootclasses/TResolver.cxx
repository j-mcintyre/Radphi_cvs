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

//////////////////////////////////////////////////////////////////  
//  This class starts w/ a list of photons from a TReconstrucion and
//  attempts to find the decay chains that lead to them.
//////////////////////////////////////////////////////////////////

#include <ieeefp.h>
#include <iostream.h>

#include <TCollection.h>
#include <TClonesArray.h>
#include <TSortedList.h>

#include "TChannel.h"
#include "TReconstruction.h"
#include "TResolver.h"

ClassImp(TResolver)

TResolver::TResolver(const TReconstruction *event)
          :fChannelLists((Int_t)kNumChannels)
{
  //constructor for a TResolver.  uses the photons
  // in the TReconstruction as its starting photons
  const Int_t numChannels=(const Int_t)kNumChannels;
  for(Int_t channelNum=0;channelNum<numChannels;channelNum++) {
    if (channelNum==(Int_t)kGamma1) {
      TClonesArray *phots= new TClonesArray("TChannel",event->GetNumPhotons());
      fChannelLists[channelNum]=phots;
      for(Int_t photNum=0;photNum<event->GetNumPhotons();photNum++) {
        const Float_t *momPhoton=event->GetPhoton(photNum);
        if (!isnanf(momPhoton[0])&&!isnanf(momPhoton[1])&&
            !isnanf(momPhoton[2])&&!isnanf(momPhoton[3]))
          new ((*phots)[photNum]) TChannel(kGamma1,momPhoton);
        if (event->GetPhoton(photNum)[0]==0.0) {
          cout << "Null photon: photNum= " << photNum << endl;
        }
      }
    }
    else 
      fChannelLists[channelNum] = new TSortedList;
  }
}

TResolver::~TResolver()
{
  // destructor.  Frees new'ed memory
  TIter iter(&fChannelLists);
  TCollection *channelList;
  while(channelList= (TCollection*)iter()) channelList->Delete();
  fChannelLists.Delete();
}

void TResolver::Display() const
{
  // prints out the channel
  TIter iter(&fChannelLists);
  TCollection *channelList;
  while (channelList = (TCollection*)iter()) {
    if (channelList->IsEmpty()) continue;
    TIter iter2(channelList);
    TChannel *channel;
    while (channel = (TChannel*)iter2()) channel->Display();
  }
}

void TResolver::FindAll(EChannelKind kind,EChannelKind from1,EChannelKind from2)
{
  // finds all kind's that can be formed from from1 & from2
  // assumes that this is a proper channel path, i.e.
  // you aren't trying to from a phi from 2 gammas.
  TCollection *channelListDest=(TCollection*)fChannelLists[kind];
  TCollection *channelListSrc1=(TCollection*)fChannelLists[from1];
  TCollection *channelListSrc2=(TCollection*)fChannelLists[from2];
  channelListDest->Delete();
  TChannel *parent;
  TIter iter1(channelListSrc1);
  TChannel *child1;
  while (child1= (TChannel*)iter1()) {
    TIter iter2(channelListSrc2);
    TChannel *child2;
    while ((child2= (TChannel*)iter2())&&(child1!=child2)) {
      parent=new TChannel(kind,*child1,*child2);
      if (parent->IsInRange()&&parent->IsConsistent()) 
        channelListDest->Add(parent);
      else
        delete parent;
    }
  }
}

void TResolver::FindAll(EChannelKind kind,EChannelKind from1,EChannelKind from2,
					  EChannelKind from3)
{
  // finds all kinds that can be formed from from1, from2, & from3
  // assumes that a kind can be formed from a from1, a from2 , and a from3,
  // i.e. you aren't making a pi0 from a phi, a eta, and an omega
  // also assumes that if from1==from3, then from2==from1 as well
  // i.e. you don't call it like FindAll(kEtaprime6,kPi02,kEta2,kPi02) ! Bad !
  // this is the right way: FindAll(kEtaprime6,kPi02,kPi02,kEta2) ! Good !
  // or you could do this: FindAll(kEta6,kPi02,kPi02,kPi02) ! Good !
  if (kind==kGamma1) 
    cout << "warning: attempt to find the decay modes of a photon" << endl;
  TCollection *channelListDest=(TCollection*)fChannelLists[kind];
  TCollection *channelListSrc1=(TCollection*)fChannelLists[from1];
  TCollection *channelListSrc2=(TCollection*)fChannelLists[from2];
  TCollection *channelListSrc3=(TCollection*)fChannelLists[from3];
  channelListDest->Delete();
  TChannel *parent;
  TIter iter1(channelListSrc1);
  TChannel *child1;
  while (child1= (TChannel*)iter1()) {
    TIter iter2(channelListSrc2);
    TChannel *child2;
    while ((child2= (TChannel*)iter2())&&(child1!=child2)) {
      if (!child2->IsCompatible(child1)) continue;
      TIter iter3(channelListSrc3);
      TChannel *child3;
      while ((child3= (TChannel*)iter3())&&(child2!=child3)) {
        parent=new TChannel(kind,*child1,*child2,*child3);
        if (parent->IsInRange()&&parent->IsConsistent()) {
          channelListDest->Add(parent);
        }
        else
          delete parent;
      }
    }
  }
}

void TResolver::FindAll(EChannelKind kind)
{
  // Finds all kind's it can, by calling the correct function
  // with the correct parameters
  switch(kind) {
    case kGamma1 : cout << "Warning: kind=kGamma1 not valid parameter \
                           in TResolver::EParticleKind" << endl;break;
    case kPi02 : FindAll(kind,kGamma1,kGamma1);break;
    case kK0L6 : FindAll(kind,kPi02,kPi02,kPi02);break;
    case kK0S4 : FindAll(kind,kPi02,kPi02);break;
    case kEta2 : FindAll(kind,kGamma1,kGamma1);break;
    case kEta6 : FindAll(kind,kPi02,kPi02,kPi02);break;
    case kOmega3Pi0 : FindAll(kind,kGamma1,kPi02);break;
    case kOmega3Eta : FindAll(kind,kGamma1,kEta2);break;
    case kOmega7 : FindAll(kind,kGamma1,kEta6);break;
    case kEtaprime2 : FindAll(kind,kGamma1,kGamma1);break;
    case kEtaprime4 : FindAll(kind,kGamma1,kOmega3Pi0);break;
    case kEtaprime6 : FindAll(kind,kPi02,kPi02,kEta2);break;
    case kf04 : FindAll(kind,kPi02,kPi02);break;
    case ka04 : FindAll(kind,kPi02,kEta2);break;
    case ka08 : FindAll(kind,kPi02,kEta6);break;
    case kPhi3Eta : FindAll(kind,kGamma1,kEta2);break;
    case kPhi3Pi0 : FindAll(kind,kGamma1,kPi02);break;
    case kPhi3Etaprime : FindAll(kind,kGamma1,kEtaprime2);break;
    case kPhi4 : FindAll(kind,kGamma1,kOmega3Pi0);break;
    case kPhi5Pi0 : FindAll(kind,kGamma1,kPi02,kPi02);break;
    case kPhi5f0 : FindAll(kind,kGamma1,kf04);break;
    case kPhi5Eta : FindAll(kind,kGamma1,kPi02,kEta2);break;
    case kPhi5a0 : FindAll(kind,kGamma1,ka04);break;
    case kPhi5Etaprime : FindAll(kind,kGamma1,kPi02,kEtaprime2);break;
    case kPhi7Eta : FindAll(kind,kGamma1,kEta6);break;
    case kPhi7Etaprime : FindAll(kind,kGamma1,kEtaprime6);break;
    case kPhi8 : FindAll(kind,kGamma1,kOmega7);break;
    default : cout << "Waring: unknown channel type in TResolver:FindAll." << endl;
              break;
  }
}

void TResolver::FindAll()
{
  //calls all the FindAll routines and finds possible channels
  FindAll(kPi02);
  FindAll(kK0L6);
  FindAll(kK0S4);
  FindAll(kEta2);
  FindAll(kEta6);
  FindAll(kOmega3Pi0);
  FindAll(kOmega3Eta);
  FindAll(kOmega7);
  FindAll(kEtaprime2);
  FindAll(kEtaprime4);
  FindAll(kEtaprime6);
  FindAll(kf04);
  FindAll(ka04);
  FindAll(ka08);
  FindAll(kPhi3Eta);
  FindAll(kPhi3Pi0);
  FindAll(kPhi3Etaprime);
  FindAll(kPhi4);
  FindAll(kPhi5Pi0);
  FindAll(kPhi5f0);
  FindAll(kPhi5Eta);
  FindAll(kPhi5a0);
  FindAll(kPhi5Etaprime);
  FindAll(kPhi7Eta);
  FindAll(kPhi7Etaprime);
  FindAll(kPhi8);
}
