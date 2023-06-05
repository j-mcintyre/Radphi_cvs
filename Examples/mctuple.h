/*
   This header file is the definition of the paw ntuple for the Radphi
   Monte Carlo simulated event. It contains the input Monte Carlo
   data as its basic block and is meant to be extended with higher-
   level analysis blocks as they are developed.  East time the format
   is changed, please update the format number identified in the title.
   This ntuple format can also hold real data, after conversion from
   the raw data to values with physical units has been completed.
   It also has version infirmation in it for compatibility reasons
*/

#include <ntypes.h>

#define NT_TITLE "Radphi MC event Ntuple version 2.11"

#define NT_FORM_HEAD  "runNo:I,eventNo:I,isMC:I"

typedef struct {
    int runNo; /* this will always be 1 */
    int eventNo;
    int isMC;
} nt_head_t;

#define NT_MAX_MC 30
#define NT_FORM_MC   "nI[0,30]:I,type:I,wgt:R,kindI(nI):I,momI(4,nI):R,\
                      nV[0,30]:I,kindV(nV):I,origV(4,nV):R,\
                      nF[0,30]:I,kindF(nF):I,vertF(nf):I,momF(4,nF):R"
typedef struct {
    int nI;
    int type;
    float wgt;
    int kindI[NT_MAX_MC];
    float momI[NT_MAX_MC][4];
    int nV;
    int kindV[NT_MAX_MC];
    float origV[NT_MAX_MC][4];
    int nF;
    int kindF[NT_MAX_MC];
    int vertF[NT_MAX_MC];
    float momF[NT_MAX_MC][4];
} nt_mc_event_t;

typedef struct    
{                      /* obsolete format for old Monte Carlo files */
  int32 channel;
  float32 energy[2];
  uint32 time[2];
}bgv_hit_old_t;

typedef struct
{
  uint32 nhits;
  bgv_hit_old_t hit[1];
}bgv_hits_old_t;

void declare_mc_ntuple(int id, nt_head_t *p_head);
void appendnt_mc_event(int id, nt_mc_event_t *p_mc);
