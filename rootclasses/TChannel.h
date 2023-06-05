
#ifndef TChannel_h
#define TChannel_h

/* the naming scheme for particles is as follows:
 *   particleName(# of photons in final state)[dominant meson]
 * the dominant meson is only listed in cases where there would 
 * otherwise be ambiguity.
 */

enum EChannelKind
{
  kGamma1,
  kPi02,
  kK0L6,kK0S4,
  kEta2,kEta6,
  kOmega3Pi0,kOmega3Eta,kOmega7,
  kEtaprime2,kEtaprime4,kEtaprime6,
  kf04,
  ka04,ka08,
  kPhi3Eta,kPhi3Pi0,kPhi3Etaprime,kPhi4,kPhi5Pi0,kPhi5f0,kPhi5Eta,
  kPhi5a0,kPhi5Etaprime,kPhi7Eta,kPhi7Etaprime,kPhi8,
  kNumChannels,kUnknownChannel
};

#include <TLorentzVector.h>
#include <TList.h>

class TObject;

class TChannel : public TObject
{

private:
  // note: TChannel uses a +--- metric, not the ROOT ---+ metric!
  TLorentzVector fp; // the momentum of the parent
  EChannelKind fKind; // the kind of decay channel
  TList fChildren; // the children of this decay
  Float_t fError; // the error of this step in the decay
  Float_t fMass; // the invariant mass of the sum of this's children's momenta
  static Float_t fRange; // the range of error that is considered plausible

  void FindError();

public:
  TChannel(const TChannel &channel);
  TChannel(EChannelKind kind,TChannel &child1,TChannel &child2);
  TChannel(EChannelKind kind,TChannel &child1,TChannel &child2,
			     TChannel &child3);
  TChannel(EChannelKind kind,TLorentzVector &p) : 
          fKind(kind),fp(p) {FindError();}
  TChannel(EChannelKind kind,const Float_t *p) :
          fKind(kind),fp(p[1],p[2],p[3],p[0]) {FindError();}
  TChannel(EChannelKind kind,Float_t E,Float_t px,Float_t py, Float_t pz):
	  fKind(kind),fp(px,py,pz,E) {FindError();}

  void SetMomentum(Float_t p[4]) {fp.SetXYZT(p[1],p[2],p[3],p[0]);FindError();}
  void SetMomentum(const TLorentzVector *p) {fp=*p;FindError();}
  void SetMomentum(Float_t E,Float_t px,Float_t py,Float_t pz)
	{fp.SetXYZT(px,py,pz,E);FindError();}

  static void SetRange(Float_t range) {fRange=range;}

  Int_t Compare(TObject *particle); // returns 1 if this has more error
				    // than particle, if = 0, if < -1
  Bool_t IsSortable() const {return kTRUE;}
  Int_t IsInRange() const {return Err()<fRange;}
  Bool_t IsCompatible(const TChannel *p) const ; //returns true if they share no children
  Bool_t IsConsistent() const; //returns true if children are compatible

  Float_t Err() const {return fError;}

  Double_t Mass() const {return M();}
  Double_t M() const {return fMass;}
  Double_t Theta() const {return(fp.Theta());}
  Double_t Phi() const {return(fp.Phi());}
  TVector3 BoostVector() const {return(fp.BoostVector());}
  void Boost(const TVector3& b) {fp.Boost(b);}

  Float_t M0() const {return M0(fKind);}  // the central mass of this particle
  Float_t Sigma() const {return Sigma(fKind);} 
				// the mass resolution of this particle
  const char *Name() const {return Name(fKind);}// get the name of this particle
  Float_t GetRange() const {return GetRange(fKind);}
  Float_t LowerBound() const {return LowerBound(fKind);}
  Float_t UpperBound() const {return UpperBound(fKind);}

  static Float_t M0(EChannelKind kind); // the central mass of kind
  static Float_t Sigma(EChannelKind kind); // the mass resolution of kind
  static const char *Name(EChannelKind kind); // get the name of kind
  static Float_t GetRange(EChannelKind kind) {return fRange*Sigma(kind);}
  static Float_t LowerBound(EChannelKind kind) {return M0(kind)-GetRange(kind);}
  static Float_t UpperBound(EChannelKind kind) {return M0(kind)+GetRange(kind);}

  void Adopt(TChannel *child) {fChildren.Add(child);}
  TChannel *GetChild(Int_t childNum) {return (TChannel*)fChildren.At(childNum);}
  TCollection *GetChildren() {return &fChildren;}
  Int_t GetNumChildren() const 
        {Int_t nc=fChildren.GetSize();return nc<0?0:nc;}

  void Display(const TString &indent="") const;

  ClassDef(TChannel,1) // decay channel class
};

#endif
