
//---Author: Tom Bogue
//---Version: 1.00/00

/*************************************************************************
 * Copyright(c) 1999, University of Connecticut, All rights reserved.    *
 * Author: Tom Bogue							 *
 *                                                                       *
 * Permission to use, copy, modify and distribute this software and its  *
 * documentation for non-commercial purposes is hereby granted without   *
 * fee, provided that the above copyright notice appears in all copies   *
 * and that both the copyright notice and this permission notice appear  *
 * in the supporting documentation. The author makes no claims about the *
 * suitability of this software for any purpose.                         *
 * It is provided "as is" without express or implied warranty.           *
 *************************************************************************/

////////////////////////////////////////////////////////////
// This class describes a particular decay channel	  //
// please note that it uses a +--- metric,		  //
// not following the ROOT std.				  //
////////////////////////////////////////////////////////////

#include "TChannel.h"
#include <iostream.h>

ClassImp(TChannel)

Float_t TChannel::fRange=2.0;

TChannel::TChannel(const TChannel &channel)
{
  //copy constructor
  fKind=channel.fKind;
  fp=channel.fp;
  fError=channel.fError;
  TIter iter(&channel.fChildren);
  TChannel *child;
  while (child= (TChannel*)iter()) fChildren.Add(child);
}

TChannel::TChannel(EChannelKind kind,TChannel &child1,TChannel &child2)
{
  // creates this as the parent of child1 and child2
  fKind=kind;
  fp=child1.fp+child2.fp;
  fChildren.Add(&child1);
  fChildren.Add(&child2);
  FindError();
}

TChannel::TChannel(EChannelKind kind,TChannel &child1,TChannel &child2,
				     TChannel &child3)
{
  // creates this as the parent of child1, child2, and child3
  TChannel(kind,child1,child2);
  fp+=child3.fp;
  fChildren.Add(&child3);
  FindError();
}

void TChannel::FindError()
{
  //finds the error, and sets the mass
  fMass=fp.M();
  Float_t m0=M0();
  fError=TMath::Abs(fMass-m0)/Sigma();
}

Int_t TChannel::Compare(TObject *channel)
{
  // compares 2 channels by error.
  // if this has less error than channel, returns -1
  // if this has just as much error as channel, returns 0
  // if this has more error than channel, returns 1
  Float_t diffErr=Err()-((TChannel *)channel)->Err();
  if (diffErr>0) return(1);
  if (diffErr<0) return(-1);
  return(0);
}

Bool_t TChannel::IsCompatible(const TChannel *p) const
{
  // returns true if this and p share no final state particles
  //  (and therefore share no particles)
  if (!fChildren.IsEmpty()) { //branch down this channel
    TIter iter(&fChildren);
    TChannel *child;
    while (child= (TChannel*)iter()) if (!child->IsCompatible(p)) return kFALSE;
  }
  else { //this is a stable state, is p?
    if (!p->fChildren.IsEmpty()) {
      TIter iter(&p->fChildren);
      TChannel *child;
      while (child= (TChannel*)iter()) if (!IsCompatible(child)) return kFALSE;
    }
    else { // OK, both this and p are stable states
      if (this==p)
        return kFALSE;
      else 
        return kTRUE;
    }
  }
  return kTRUE;
}

Bool_t TChannel::IsConsistent() const
{
  //returns true if this's children are compatible
  TIter iter1(&fChildren);
  TChannel *child1;
  while (child1= (TChannel*)iter1()) {
    TIter iter2(&fChildren);
    TChannel *child2;
    while ((child2= (TChannel*)iter2())&&(child2!=child1)) {
      if (!child1->IsCompatible(child2)) return kFALSE;
    }
  }
  return kTRUE;
}

Float_t TChannel::M0(EChannelKind kind)
{
  // returns the true mass that this particle should have
  switch(kind) {
    case kGamma1 : return 0.0 ;break;
    case kPi02 : return 0.1349764 ;break;
    case kK0L6 : 
    case kK0S4 : return 0.497672;break;
    case kEta2 : 
    case kEta6 : return 0.5473;break;
    case kOmega3Pi0 :
    case kOmega3Eta :
    case kOmega7 : return 0.78194 ;break;
    case kEtaprime2 :
    case kEtaprime4 :
    case kEtaprime6 : return 0.95778 ;break;
    case kf04 : return 0.980;break;
    case ka04 :
    case ka08 : return 0.9834;break;
    case kPhi3Eta :
    case kPhi3Pi0 :
    case kPhi3Etaprime :
    case kPhi4 :
    case kPhi5Pi0 :
    case kPhi5f0 :
    case kPhi5Eta :
    case kPhi5a0 :
    case kPhi5Etaprime :
    case kPhi7Eta :
    case kPhi7Etaprime :
    case kPhi8 : return 1.019413;break;
    default : cout << "Warning: unknown particle in TChannel::Mass. kind= " << (Int_t) kind << endl;
              return -1.0 ;break;
  }
}

Float_t TChannel::Sigma(EChannelKind kind)
{
  // returns the experimental spread in mass for this particle
  switch(kind) {
    case kGamma1 : return 0.01;break;
    case kPi02 : return 0.0150464;break;
    case kK0L6 : 
    case kK0S4 : return 0.03;break;
    case kEta2 : return 0.041;break;
    case kEta6 :
    case kOmega3Pi0 :
    case kOmega3Eta :
    case kOmega7 :
    case kEtaprime2 :
    case kEtaprime4 :
    case kEtaprime6 :
    case kf04 :
    case ka04 :
    case ka08 :
    case kPhi3Eta :
    case kPhi3Pi0 :
    case kPhi3Etaprime :
    case kPhi4 :
    case kPhi5Pi0 :
    case kPhi5f0 :
    case kPhi5Eta :
    case kPhi5a0 :
    case kPhi5Etaprime :
    case kPhi7Eta :
    case kPhi7Etaprime :
    case kPhi8 : return 0.06;
    default : cout << "Warning: unknown particle in TChannel::Sigma kind " << endl;
              return -1.0 ;break;
  }
}

const char *TChannel::Name(EChannelKind kind)
{
  // returns the name of this decay channel
  switch(kind) {
    case kGamma1 :	return "gamma1";break;
    case kPi02 :	return "pi0";break;
    case kK0L6 :	return "K0L6";break;
    case kK0S4 :	return "K0S4";break;
    case kEta2 :	return "eta2";break;
    case kEta6 :	return "eta6";break;
    case kOmega3Pi0 :	return "omega3Pi0";break;
    case kOmega3Eta :	return "omega3Eta";break;
    case kOmega7 :	return "omega7";break;
    case kEtaprime2 :	return "Eta'2";break;
    case kEtaprime4 :	return "Eta'4";break;
    case kEtaprime6 :	return "Eta'6";break;
    case kf04 :		return "f04";break;
    case ka04 :		return "a04";break;
    case ka08 :		return "a08";break;
    case kPhi3Eta :	return "phi3Eta";break;
    case kPhi3Pi0 :	return "phi3Pi0";break;
    case kPhi3Etaprime :return "phi3Eta'";break;
    case kPhi4 :	return "phi4";break;
    case kPhi5Pi0 :	return "phi5Pi0";break;
    case kPhi5f0 :	return "phi5f0";break;
    case kPhi5Eta :	return "phi5Eta";break;
    case kPhi5a0 :	return "phi5a0";break;
    case kPhi5Etaprime :return "phi5Eta'";break;
    case kPhi7Eta :	return "phi7Eta";break;
    case kPhi7Etaprime :return "phi7Eta'";break;
    case kPhi8 :	return "phi8";break;
    default : cout << "Warning: unknown channel in TChannel::Name" << endl;
              return "" ;break;
  }
}

void TChannel::Display(const TString &indent) const
{
  // prints out this decay channel
  cout << indent << Name() << " w/ p= (" << fp.E() << ',' << fp.Px() 
       << ',' << fp.Py() << ',' << fp.Pz() << ")"  
       << " Err= " << Err() << endl;
  if (GetNumChildren()>0) {
    TString newIndent=indent+"  ";
    TIter iter(&fChildren);
    TChannel *child;
    while (child= (TChannel*)iter()) child->Display(newIndent);
  }
}
