/*
   This header file is the definition of the paw ntuple for the
   Radphi scalers.  It contains the scaler data as the basic block
   and is meant to be extended with additional scalers as needed.
   Each time the format is changed, please update the format number
   identified in the title.
*/

#define SC_TITLE "Radphi scaler Ntuple version 2.00"

#define FORM_HEAD "run:I,event:I,seq:I"

#define BUFSIZE 100000
#define NTID 800
#define MAX_SCALERS 1024

typedef struct scaler_label {
  const char *string;
  const char *tag;
} scalabel;

typedef struct scaler_table {
  int nscalers;
  scalabel *scaler[MAX_SCALERS];
} scatable;

void declare_scaler_ntuple (int id, scatable *names, int *p_head, int *p_scalers);
