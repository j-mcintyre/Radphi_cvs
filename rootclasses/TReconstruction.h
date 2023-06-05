#ifndef TReconstruction_h
#define TReconstruction_h

#include "TChannel.h"

class TTree;
class TH1;
class TH2;

class TReconstruction
{
private :
   TTree          *fTree; //pointer to the analyzed TTree
//Declaration of leaves types
   Int_t           nphot; // the total # of photons found
   Int_t	   nfrwd; // the # of photons found in the forward detector(the LGD)
   Float_t         pvect[400][4]; // the momentum vectors of these photons

public :
   TReconstruction(TTree *tree=0);
   ~TReconstruction() {;}
   Int_t GetEvent(Int_t event);
   void Init(TTree *tree);
   void Show(Int_t event = -1) const;
   Int_t GetNumPhotons() const {return(nphot);}
   const Float_t *GetPhoton(Int_t photonNum) const {return(pvect[photonNum]);}
   void FindPi0(Int_t eventNum);
   void FindPi0();
   TH1 *FitParticle(EChannelKind particle);
   TH1 *PlotMof2();
   TH2 *Plot2Particles();
   TH2 *PlotThetaPhiOfPi0();
   TH1 **PlotEBGV_ELGD();
   TH1 *PlotChannelMass(EChannelKind kind);
   TH1 **FindSigma(Int_t maxEntries=999999999);
   TH1 *PlotMofBGVphoton();
   TH1 *PlotThetaOfGammas();

   ClassDef(TReconstruction,1) // Reconstruction data class
};

#endif
