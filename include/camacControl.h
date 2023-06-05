/*
  $Log: camacControl.h,v $
  Revision 1.6  2000/05/02 15:10:45  radphi
  Changes made by radphi@urs3
  adding multi-crate, and turning off controller inhibit line

  Revision 1.5  1999/06/22 03:16:31  radphi
  Changes made by radphi@urs1
  slowly but surely updating the function prototypes so that the init function
  actually takes a crate number, rather than it being hardwired to one.

  Revision 1.4  1997/05/17 21:43:23  radphi
  Changes made by lfcrob@jlabs2
  Added camac_disconnect()

 * Revision 1.3  1997/05/16  01:41:09  radphi
 * Changes made by lfcrob@jlabs2
 * Added lrs 2552 routines
 *
 * Revision 1.2  1997/05/08  05:02:16  radphi
 * Changes made by lfcrob@dustbunny
 * Added energy sum
 *
 * Revision 1.1  1997/04/08  21:07:39  radphi
 * Initial revision
 *
  */

#ifndef CAMAC_CONT_H_INCLUDED
#define CAMAC_CONT_H_INCLUDED

#define CAMAC_OK 0
#define CAMAC_ERROR 1
#define CAMAC_CAMAC_ERROR 2

/* For triggerBase.c */

int camac_init(void);
int camac_initStatus(void);
int camac_getModule(char *alias,int *myCrate, int *mySlot);
int camac_setInhibit(int inhibit);
void camac_disconnect(void);


int lrs3420_init(int crate, int slot);
int lrs3420_writeThreshold(int slot, int channel, int theshold);
int lrs3420_writeThresholds(int slot, int threshold[16]);
int lrs3420_readThreshold(int slot, int channel, int *theshold);
int lrs3420_readThresholds(int slot, int threshold[16]);
int lrs3420_writeMask(int slot, int mask);
int lrs3420_readMask(int slot, int *mask);
int lrs3420_writeDuration(int slot, int duration, int deadtime);
int lrs3420_readDuration(int slot, int *duration, int *deadtime);
int lrs3420_testPulse(int slot);
int doLRS_3420(int crate, int slot);


int lrs4413_init(int crate, int slot);
int lrs4413_writeThreshold(int slot, int threshold);
int lrs4413_setToManualThreshold(int slot);
int lrs4413_readThreshold(int slot, int *threshold);
int lrs4413_writeMask(int slot, int mask);
int lrs4413_readMask(int slot, int *mask);
int lrs4413_localMode(int slot);
int lrs4413_remoteMode(int slot);
int doLRS_4413(int crate, int slot);

int lrs4516_init(int slot);
int lrs4516_setOR(int slot,int stage);
int lrs4516_setAND(int slot,int stage);
int lrs4516_test(int slot,int stage);
int doLRS_4516(int crate, int slot);

int iuVarDelay_init(int slot);
int iuVarDelay_write(int slot, int channel,int delay);
int iuVarDelay_read(int slot, int channel,int *delay);
int doIU_delay(int crate, int slot);

int iuESum_init(int slot);
int iuESum_write(int slot,int threshold);
int iuESum_read(int slot,int *threshold);
int doIU_ESum(int crate, int slot);

int lrs2323_init(int slot);
int lrs2323_setGate(int slot, int channel, int mantissa, int exponent,int latchMode, int delayWidth);
int lrs2323_readGate(int slot, int channel, int *mantissa, int *exponent,int *latchMode, int *delayWidth);
int lrs2323_stopGate(int slot, int channel);
int lrs2323_startGate(int slot, int channel);
int doLRS_2323(int crate, int slot);

int jorway41_init(int slot);
int jorway41_writeLevels(int slot, int pattern);
int jorway41_writePulses(int slot, int pattern);
int jorway41_selectWrite(int slot, int pattern);
int jorway41_selectClear(int slot, int pattern);
int doJorway_41(int crate,int slot);

int ks3922_init(int crate);
int ks3922_remove_inhibit(int crate);
int doKS_3922(int crate);

int lgdXP_init(void);
int lgdXP_status(int *status);
int lgdXP_startMove(int motor);
int lgdXP_stopMove(int motor);
float lgdXP_getIPS(int direction);

int joergerSMCR_init(int crate, int slot);
int joergerSMCR_status(int slot, int *status);
int joergerSMCR_move(int slot,int motor,int direction,int nSteps);
 

#define JOERGER_SMCR_CW 0x1
#define JOERGER_SMCR_CCW 0x2

#define JOERGER_SMCR_MAX_STEPS 0x7fff
#define JOERGER_SMCR_CCW_MASK 0x8000


#define JOERGER_SMCR_M0_CW_LIMIT 0x1
#define JOERGER_SMCR_M0_CCW_LIMIT 0x2
#define JOERGER_SMCR_M0_DONE 0x4
#define JOERGER_SMCR_M1_CW_LIMIT 0x8
#define JOERGER_SMCR_M1_CCW_LIMIT 0x10
#define JOERGER_SMCR_M1_DONE 0x20
#define JOERGER_SMCR_NO_POWER 0x40
#define JOERGER_SMCR_POWER_STUCK 0x80

#define KS_3922_CONTROL_SLOT 30

int ks3516_init(int crate, int slot);
int ks3516_writeDMA(int slot,int channel);
int ks3516_writeCMA(int slot,int channel);
int ks3516_readData(int slot,int channel, int *value);
int ks3516_readGain(int slot,int channel, int *value);
int ks3516_writeGain(int slot,int channel, int gain);
int ks3516_startScan(int slot);
int ks3516_stopScan(int slot);
int doKS_3516(int crate,int slot);


int lgdHV_init(void);
int lgdHV_ready(void);
int lgdHV_setChannel(int row, int col, int voltage);
int lgdHV_setRow(int row, int voltage);
int lgdHV_setCol(int col, int voltage);
int lgdHV_setAll(int voltage);
int lgdHV_pulseChannel(int row, int col);
int lgdHV_pulseRow(int row);
int lgdHV_pulseCol(int col);
int lgdHV_pulseAll(void);

#define LGDHV_READY 0
#define LGDHV_BUSY 3

int lrs2551_init(int slot);
int lrs2551_readChannel(int slot, int channel, int *data);
int lrs2551_clear(int slot);
int lrs2551_increment(int slot);
int lrs2551_readAll(int slot, int data[12]);

#define L3_E_OR_M 1
#define L3_E_AND_M 0
int lrs2372_init(int crate, int slot);
int lrs2372_writeMode(int slot, int mode);
int lrs2372_writeAddress(int slot, int address);
int lrs2372_writeData(int slot, int L3_decision_operator, 
		      int E_thresh, int M_thresh);

#endif
