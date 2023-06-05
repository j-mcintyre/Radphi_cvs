#ifndef TMCdata_h
#define TMCdata_h

#include <TTree.h>
#include <TFile.h>

class TMCdata
{
private :
   TTree          *fTree; //pointer to the analyzed TTree
//Declaration of leaves types
   Int_t           runno; // the run number of this event
   Int_t           eventno; //the event number of this event
   Int_t           ismc; // a flag to say if there is MC data in this run

   Int_t           ni; // the number of initial state particles
   Int_t           type; // a user-defined type of this channel
   Float_t         wgt; // the statistical weight for this event
   Int_t           kindi[30]; // the kind of initial state particle[i]
   Float_t         momi[30][4]; // the momentum 4-vector of particle[i]
   Int_t           nv; // the number of final state vertices
   Int_t           kindv[30]; // the kind of final state vertex[i]
   Float_t         origv[30][4]; // the (t,x,y,z) origin of vertex[i]
   Int_t           nf; // the # of final state particles
   Int_t           kindf[30]; // the kind of final state particle[i]
   Int_t           vertf[30]; // vertex index of the final particle[i]
   Float_t         momf[30][4]; // the momentum of the final particle[i]

public:
   TMCdata(TTree *tree=0); 
   ~TMCdata() {;}
   Int_t GetEvent(Int_t event);
   void Init(TTree *tree);
   void Show(Int_t event = -1) const;
   
   ClassDef(TMCdata,1) // MC data class
};

#endif
