/*
   This header file is the definition of the paw ntuple for the
   raw Radphi event.  It contains the raw data as the basic block
   and is meant to be extended with higher-level analysis blocks
   as they are developed.  Each time the format is changed,
   please update the format number identified in the title.
*/

#define NT_TITLE "Radphi event Ntuple version 2.12"
#define NT_FORM_HEAD "runNo:I, eventNo:I, isMC:I"

#define NT_FORM_RAW "evType:I, meMAM:I, cMAM:I, eMAM:I, m2MAM:I,\
                     nBSD[0,400]:I, iBSD(nBSD):I, aBSD(nBSD):I, tBSD(nBSD):I,\
                     nBGV[0,400]:I, iBGV(nBGV):I, aBGV(nBGV):I, tBGV(nBGV):I,\
                     nCPV[0,250]:I, iCPV(nCPV):I, aCPV(nCPV):I, tCPV(nCPV):I,\
                     nUPV[0,100]:I, iUPV(nUPV):I, aUPV(nUPV):I, tUPV(nUPV):I,\
                     nLGD[0,640]:I, iLGD(nLGD):I, jLGD(nLGD):I, aLGD(nLGD):I,\
                     nTAG[0,320]:I, iTAG(nTAG):I, tTAG(nTAG):I"
   
#define NT_MAX_BSD_RAW 400
#define NT_MAX_BGV_RAW 400
#define NT_MAX_CPV_RAW 250
#define NT_MAX_UPV_RAW 100
#define NT_MAX_LGD_RAW 640
#define NT_MAX_TAG_RAW 320

typedef struct {
    int runNo;
    int eventNo;
    int isMC;
} nt_head_t;

typedef struct ntuple_raw {
    int evType;
    int meMAM;
    int cMAM;
    int eMAM;
    int m2MAM;
    int nBSD;
    int iBSD[NT_MAX_BSD_RAW];
    int aBSD[NT_MAX_BSD_RAW];
    int tBSD[NT_MAX_BSD_RAW];
    int nBGV;
    int iBGV[NT_MAX_BGV_RAW];
    int aBGV[NT_MAX_BGV_RAW];
    int tBGV[NT_MAX_BGV_RAW];
    int nCPV;
    int iCPV[NT_MAX_CPV_RAW];
    int aCPV[NT_MAX_CPV_RAW];
    int tCPV[NT_MAX_CPV_RAW];
    int nUPV;
    int iUPV[NT_MAX_UPV_RAW];
    int aUPV[NT_MAX_UPV_RAW];
    int tUPV[NT_MAX_UPV_RAW];
    int nLGD;
    int iLGD[NT_MAX_LGD_RAW];
    int jLGD[NT_MAX_LGD_RAW];
    int aLGD[NT_MAX_LGD_RAW];
    int nTAG;
    int iTAG[NT_MAX_TAG_RAW];
    int tTAG[NT_MAX_TAG_RAW];
} nt_raw_event_t;

void declare_raw_ntuple(int id, nt_head_t *p_head);
void appendnt_mc_event(int id, nt_raw_event_t *p_raw);
void appendnt_raw_event(int id, nt_raw_event_t *p_raw);
