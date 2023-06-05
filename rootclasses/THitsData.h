#ifndef THitsData_h
#define THitsData_h

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>

class THitsData
{
private :
   TTree          *fTree; //pointer to the analyzed TTree
//Declaration of leaves types
   Int_t           runno; // the run number of this event
   Int_t           eventno; // the event number of this event
   Int_t           ismc; // a flag to state if there is MC data for this event

   Int_t           nhbsd; // the number of hits in the BSD for this event
   Int_t           chbsd[50]; // the channels of those hits
   Float_t         Ebsd[50]; // the energy of those hits
   Int_t           ntbsd[50]; // the number of TDC's for each hit
   Int_t           t1bsd[50]; // the index in the timelist of the 1st TDC for this hit

   Int_t           nhbgv; // the number of hits in the BGV for this event
   Int_t           chbgv[50]; // the channel of those hits
   Float_t         Ebgvdwn[50]; // the energy at the downstream ADC of each hit
   Float_t         Ebgvup[50]; // the energy at the upstream ADC of each hit
   Float_t         tbgvdwn[50]; // the downsteam TDC for each hit
   Float_t         tbgvup[50]; // the upstream TDC for each hit

   Int_t           nhcpv; // the number of hits in the CPV for this event
   Int_t           chcpv[50]; // the channel #'s for these hits
   Float_t         Ecpv[50]; // the energy of each hit
   Int_t           ntcpv[50]; // the # of TDC's for each hit
   Int_t           t1cpv[50]; // the index in the timelist of the 1st TDC for this hit

   Int_t           nhupv; // the number of hits in the UPV for this event
   Int_t           chupv[50]; // the channel #'s of those hits
   Float_t         Eupv[50]; // the energy of each hit
   Int_t           ntupv[50]; // the # of TDC's for each hit
   Int_t           t1upv[50]; // the index in the timelist of the 1st TDC for this hit

   Int_t           nhlgd; // the number of hits in the LGD for this event
   Int_t           chlgd[500]; // the channel #'s for each hit
   Float_t         Elgd[500]; // the energy of each hit

   Int_t           ntimes; // the number of TDC's hits total for this event
   Float_t         le[500]; // the timelist, containing all TDC's for this event

public:
   THitsData(TTree *tree=0);
   ~THitsData() {;}
   Int_t GetEvent(Int_t event);
   void Init(TTree *tree);
   void Show(Int_t event = -1) const;
   TH1 *Plotnhbgv();
   
   ClassDef(THitsData,1) // hits data storage class
};

#endif
