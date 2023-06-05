/*
  $Log: unpackData.c,v $
  Revision 1.23  2003/09/05 15:46:31  radphi
  -unpackData.c - modified to do nothing and return 0 if there is no raw data
       group in the event; makes it easier to analyse both real data and MC
       with the same analysis code [rtj]

  Revision 1.22  2002/08/27 00:56:55  radphi
  - took some test code out of unpackData.c -rtj-

  Revision 1.21  2000/07/13 13:24:18  radphi
  Changes made by jonesrt@grendl
  introduced a more compact form of the raw data group for efficiency -rtj

  Revision 1.20  2000/05/29 20:42:22  radphi
  Changes made by radphi@urs3
  unpacking of MAM and DOR (upper and lower disrims) for IU ADC
  dsa

  Revision 1.19  2000/01/18 04:03:03  radphi
  Changes made by jonesrt@golem
  fixed endian-conversion package to allow different structures with the
  same group number but different subgroups -- to support MonteCarlo group. -rtj

  Revision 1.18  1999/06/17 17:33:02  radphi
  Changes made by radphi@jlabs1
   include BSD BGV information
     D.S. Armstrong

  Revision 1.17  1998/06/26 23:57:01  radphi
  Changes made by radphi@jlabs1
  added LRS 1875a decoding  D.S. Armstrong

  Revision 1.16  1998/06/22 21:16:06  radphi
  Changes made by jonesrt@jlabs4
  Standardized the lengths of adc,tdc groups to hold multi-hit tdc counts -RTJ

  Revision 1.15  1998/06/10 17:44:39  radphi
  Changes made by radphi@urs2
  reset reg_vals group   D.S. Armstrong

  Revision 1.14  1998/04/12 02:59:50  radphi
  Changes made by radphi@jlabs3
  add group creation for MAM   D.S. Armstrong

 * Revision 1.13  1998/04/02  20:04:47  radphi
 * Changes made by jonesrt@jlabs4
 * a dangling if was fixed in unpackData.c - R.T.Jones
 *
 * Revision 1.12  1998/04/02  17:03:04  radphi
 * Changes made by jonesrt@jlabs4
 * fixing the LGD readout - R.T.Jones
 *
 * Revision 1.11  1998/04/01  20:29:47  radphi
 * Changes made by jonesrt@jlabs3
 * correct a bug for events with missing 10C6 header word - R.T.Jones
 *
 * Revision 1.10  1998/04/01  17:37:25  radphi
 * Changes made by jonesrt@jlabs3
 * bringing the code up to date with the database - R.T.Jones
 *
 * Revision 1.9  1998/03/31  21:23:26  radphi
 * Changes made by radphi@jlabs2
 * added slot to raw structures
 *
 * Revision 1.8  1998/03/31  03:37:24  radphi
 * Changes made by jonesrt@zeus
 * getting the decoding working - R.T.Jones
 *
 * Revision 1.7  1998/03/31  02:59:14  radphi
 * Changes made by jonesrt@zeus
 * small bug fix - R.T.Jones
 *
 * Revision 1.6  1998/03/30  21:36:40  radphi
 * Changes made by radphi@urs3
 * update iu adc unpacking, add upv
 *
 * Revision 1.5  1998/03/30  19:36:31  radphi
 * Changes made by jonesrt@
 * new implementation for March 1998 run period - R.T.Jones
 *
 * Revision 1.4  1997/05/29  20:34:12  radphi
 * Changes made by lfcrob@dustbunny
 * Added <string.h> so it will compile clean
 *
 * Revision 1.3  1997/05/20  09:06:12  radphi
 * Changes made by lfcrob@jlabs2
 * Put in changes to properly handle scalers, ESUM, and CPV tdc
 *
 * Revision 1.2  1997/05/20  04:36:59  radphi
 * Changes made by radphi@jlabs2
 * Added hack to bypass CPV TDCS and ESUM ADC's
 *
 * Revision 1.1  1997/05/16  05:40:31  radphi
 * Initial revision by lfcrob@jlabs2
 * Master raw data unpacking routines - from rodd
 *
  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <disData.h>
#include <eventType.h>
#include <triggerType.h>
#include <libdecode.h>
#include <umap.h>
#include <unpackData.h>

static int verbose=0;
static FILE *logfp = NULL;


static adc_values_t *adc_vals=NULL;
static adc_values_t *rpd_adcs=NULL;
static adc_values_t *lgd_adcs=NULL;
static adc_values_t *cpv_adcs=NULL;
static adc_values_t *upv_adcs=NULL;
static adc_values_t *bsd_adcs=NULL;
static adc_values_t *bgv_adcs=NULL;
static adc_values_t *esum_adcs=NULL;
static tdc_values_t *tdc_vals=NULL;
static tdc_values_t *rpd_tdcs=NULL;
static tdc_values_t *cpv_tdcs=NULL;
static tdc_values_t *upv_tdcs=NULL;
static tdc_values_t *bsd_tdcs_long=NULL;
static tdc_values_t *bsd_tdcs_short=NULL;
static tdc_values_t *bgv_tdcs=NULL;
static tdc_values_t *tagger_tdcs=NULL;
static scaler_values_t *scaler_vals=NULL;
static scaler_values_t *scalers=NULL;
static reg_values_t *reg_vals=NULL;

static void init_groups(void);
static void reset_groups(void);
static void malloc_error( char* source, int line );
static void umap_error(int ROC, int slot, int channel) ;
static int dbmap(int ROCnum, adc_values_t* adc_group, 
		 tdc_values_t *tdc_group, scaler_values_t* scaler_group);

void unpackSetLogFile(FILE *fp)
{
  logfp = fp;
}
void unpackSetVerbose(int level)
{
  verbose=level;
}

int unpackEvent(itape_header_t *event, int BUFSIZE)

{
  int *rawData=NULL; /* pointer to current location in raw data */
  int rawSize; /* Total event size */
  int length=0; /* Length of current bank */
  int ERR=0;
  int ROCnum;
  int size;
  void *group;
  int typeSave;
  int ii;

  if(event->eventType != EV_DATA){
    return(0);
  }
  if((rawData = data_getGroup(event,GROUP_RAW,0)) == NULL){
    itape_header_t *frozenData = data_getGroup(event,GROUP_FROZEN,0);
    if(frozenData == NULL){
  /*
   * In the past, unpackEvent() would return with success without doing
   * anything if the event was not of EV_DATA type, but if an EV_DATA
   * event had no raw group (eg. Monte Carlo) then it would print an
   * error message and return failure.  The new semantics are for
   * unpackEvent() to assume that, if the raw group does not exist, its
   * task has already been completed and return immediately with success.
   * This allows Monte Carlo data to be processed by the same code as
   * real data, without recompiling or resorting to "if (isMC)" blocks.
      fprintf(stderr,"unpackEvent: no GROUP_RAW\n");
      return(1);
   */
      return 0;
    }
    else {
      return thawData(event, BUFSIZE);
    }
  }
  init_groups();
  rawSize = *rawData;

  /* Skip over the Event header (2 words0 and the EB bank (5 words). 
     Make the pointer point to the beginning of the first ROC bank */
  
  rawData += 7;
  rawSize -= 7;
  
  /*
    length = *rawData;
    ROCnum = ((rawData[8] >> 16) & (0xf));
    rawData = &rawData[9]; length -= 9;
  */
  while(rawSize > 0){
    length = *rawData;
    rawData++;rawSize--;
    rawSize -= length;
    ROCnum = ((*rawData) >> 16) & 0xff;
    length--;
    rawData++;
/*    fprintf(stderr,"ROC Bank : length %d, ROCnum %d , firstWord 0x%x\n",length,ROCnum,*rawData);  */
    while (length > 0) {
      switch (*rawData) {
      case START_1885:
      case START_1885_OLD:
	rawData++; length--;
	decode_1885(&rawData, &length, adc_vals);
	if ((*rawData != END_1885) && (*rawData != END_1885_OLD)) {
	  fprintf(stderr, "Didn't find end marker for 1885 %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: Did not find end marker for 1885_ADC\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n", event->runNo,
		    event->eventNo, length);
	  }
	  length = -1; ERR = 1;
	  continue;    /* takes us to top of while length loop, which is now false */
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (1885)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp){
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  }
	  length = -1; ERR = 1;
	  continue;
	}
	if(verbose){
	  if(logfp)
	    fprintf(logfp, "  Read %d channels of 1885_ADC\n", adc_vals->nadc);
	}
	rawData++; length--;
	if (dbmap(ROCnum, adc_vals, NULL, NULL)) {
	  ERR = 1; length = -1;
	  continue;
	}
	break;
	
      case START_10C6 :
	rawData++; length--;
	decode_10C6(&rawData, &length, tdc_vals);
	if (*rawData != END_10C6) {
	  fprintf(stderr, "Didn't find end marker! %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: 10C6 TDC module did not contain end marker\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n",
		    event->runNo, event->eventNo, length);
	  }
	  ERR = 1; length = -1;
	  continue;
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (10C6)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp)
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  
	  ERR = 1; length = -1;
	  continue;
	}
	if(verbose){
	  if(logfp)
	    fprintf(logfp, "  Read %d channels of 10C6_TDC\n", tdc_vals->ntdc);
	}
	rawData++; length--;
	if (dbmap(ROCnum, NULL, tdc_vals, NULL)) {
	  ERR = 1; length = -1;
	  continue;
	}
	break;

      case START_NIM_7200:
      case START_ECL_7200:
	typeSave = *rawData;
	rawData++; length--;
	decode_7200(&rawData, &length, scaler_vals);
	switch(typeSave){
	case START_NIM_7200:
	  if(*rawData != END_NIM_7200){
	    fprintf(stderr, "Didn't find end marker! %#x\n", *rawData);
	    if(logfp){
	      fprintf(logfp, "ERROR: 7200 NIM scaler module did not contain end marker\n");
	      fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n",
		      event->runNo, event->eventNo, length);
	    }
	    ERR = 1; length = -1;
	    continue;
	  }
          {  /* this scaler was hard-wired in the database to slot 5 */
            int i;
            for(i=0;i<scaler_vals->nscalers;i++)
              scaler_vals->scaler[i].slot = 5;
          }
	  break;
	case START_ECL_7200:
	  if(*rawData != END_ECL_7200){
	    fprintf(stderr, "Didn't find end marker! %#x\n", *rawData);
	    if(logfp){
	      fprintf(logfp, "ERROR: 7200 ECL scaler module did not contain end marker\n");
	      fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n",
		      event->runNo, event->eventNo, length);
	    }
	    ERR = 1; length = -1;
	    continue;
	  }
          {  /* this scaler was hard-wired in the database to slot 3 */
            int i;
            for(i=0;i<scaler_vals->nscalers;i++)
              scaler_vals->scaler[i].slot = 3;
          }
	  break;
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (7200)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp)
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  ERR = 1; length = -1;
	  continue;
	}
	if(verbose){
	  if(logfp)
	    fprintf(logfp, "  Read %d channels of 7200_scaler\n", scaler_vals->nscalers);
	}
	rawData++;length--;
	if (dbmap(ROCnum, NULL, NULL, scaler_vals)) {
	  ERR = 1; length = -1;
	  continue;
	}
	break;

      case START_1151 :
	rawData++; length--;
	decode_1151(&rawData, &length, scaler_vals);
	if (*rawData != END_1151) {
	  fprintf(stderr, "Didn't find end marker! %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: 1151 scaler module did not contain end marker\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n",
		    event->runNo, event->eventNo, length);
	  }
	  ERR = 1; length = -1;
	  continue;
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (1151)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp)
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  ERR = 1; length = -1;
	  continue;
	}
	if(verbose){
	  if(logfp)
	    fprintf(logfp, "  Read %d channels of 1151_scaler\n", scaler_vals->nscalers);
	}
	rawData++; length--;
	if (dbmap(ROCnum, NULL, NULL, scaler_vals)) {
	  ERR = 1; length = -1;
	  continue;
	}
	break;

      case START_BLOCK_7200:
	rawData++; length--;
	decode_block7200(&rawData, &length, scaler_vals);
	if (*rawData != END_BLOCK_7200) {
	  fprintf(stderr, "Didn't find end marker! %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: 7200 scaler block did not contain end marker\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n",
		    event->runNo, event->eventNo, length);
	  }
	  ERR = 1; length = -1;
	  continue;
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (7200)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp)
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  ERR = 1; length = -1;
	  continue;
	}
	if(verbose){
	  if(logfp)
	    fprintf(logfp, "  Read %d channels of 7200_scaler\n", scaler_vals->nscalers);
	}
	rawData++; length--;
	if (dbmap(ROCnum, NULL, NULL, scaler_vals)) {
	  ERR = 1; length = -1;
	  continue;
	}
	break;
	
      case START_IU_ADC:
	rawData++; length--;
 	/* in blockread IUadc data there is a special word to ignore */
	if (*rawData == ADC_BOUNDARY) {
	  rawData++; length--;
	}
	decode_IUadc(&rawData, &length, adc_vals);
	if (*rawData != END_IU_ADC) {
	  fprintf(stderr, "Didn't find end marker for IU_ADC %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: Did not find end marker for IU_ADC\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n", event->runNo,
		    event->eventNo, length);
	  }
	  length = -1; ERR = 1;
	  continue;    /* takes us to top of while length loop, which is now false */
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (IU_ADC)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp){
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  }
	  length = -1; ERR = 1;
	  continue;
	}
	if(verbose){
	  if(logfp)
	    fprintf(logfp, "  Read %d channels of IU_ADC\n", adc_vals->nadc);
	}
	rawData++; length--;
	if (dbmap(ROCnum, adc_vals, NULL, NULL)) {
	  ERR = 1; length = -1;
	  continue;
	}

	break;

      case START_1877:
	rawData++; length--;
	if(decode_1877(&rawData, &length, tdc_vals)==1){
	  fprintf(stderr, "L2 reset failure for 1877 - will be no TDC data\n");
	  if(logfp){
	    fprintf(logfp, "L2 reset failure for 1877 - will be no TDC data\n");
	  }
	}
        if (length==0) break;
	if (*rawData != END_1877) {
	  fprintf(stderr, "Didn't find end marker for 1877 %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: Did not find end marker for 1877_TDC\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n", event->runNo,
		    event->eventNo, length);
	  }
	  length = -1; ERR = 1;
	  continue;    /* takes us to top of while length loop, which is now false */
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (1877)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp){
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  }
	  length = -1; ERR = 1;
	  continue;
	}
	if(verbose){
	  if(logfp) {
	    fprintf(logfp, "  Read %d channels of 1877_TDC\n", tdc_vals->ntdc);
	    for (ii=0; ii < tdc_vals->ntdc; ii++){
	      fprintf(logfp, "      slot %d , channel %d , value  %d\n", 
   tdc_vals->tdc[ii].slot, tdc_vals->tdc[ii].channel, tdc_vals->tdc[ii].le);
	    }
	  }
	}
	rawData++; length--;
	if (dbmap(ROCnum, NULL, tdc_vals, NULL)) {
	  /*	  ERR = 1; length = -1; */
	  	  ERR = 1; 
	  continue;
	}
	break;

      case START_1875:
	rawData++; length--;
	decode_1875(&rawData, &length, tdc_vals);
	if (*rawData != END_1875) {
	  fprintf(stderr, "Didn't find end marker for 1875 %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: Did not find end marker for 1875_TDC\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n", event->runNo,
		    event->eventNo, length);
	  }
	  length = -1; ERR = 1;
	  continue;    /* takes us to top of while length loop, which is now false */
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (1875)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp){
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  }
	  length = -1; ERR = 1;
	  continue;
	}
	if(verbose){
	  if(logfp){
	    fprintf(logfp, "  Read %d channels of 1875_TDC\n", tdc_vals->ntdc);
	      for (ii=0; ii < tdc_vals->ntdc; ii++){
	      fprintf(logfp, "      slot %d , channel %d , value  %d\n", 
   tdc_vals->tdc[ii].slot, tdc_vals->tdc[ii].channel, tdc_vals->tdc[ii].le);
	    }
	  }
	}
	rawData++; length--;
	if (dbmap(ROCnum, NULL, tdc_vals, NULL)) {
	  ERR = 1; length = -1;
	  continue;
	}
	break;

      case START_IU_MAM:
	rawData++; length--;
	decode_IUmam(&rawData, &length, reg_vals);
	if (*rawData != END_IU_MAM) {
	  fprintf(stderr, "Didn't find end marker for IU_MAM %#x\n", *rawData);
	  if(logfp){
	    fprintf(logfp, "ERROR: Did not find end marker for IU_MAM\n");
	    fprintf(logfp, "       Run: %d, Event: %d, Remaining words: %d\n", event->runNo,
		    event->eventNo, length);
	  }
	  length = -1; ERR = 1;
	  continue;    /* takes us to top of while length loop, which is now false */
	}
	if (length < 0) {
	  fprintf(stderr, "Read beyond bounds of current event (IU_MAM)!\n");
	  fprintf(stderr, "Run: %d, Event: %d\n", event->runNo, event->eventNo);
	  if(logfp){
	    fprintf(logfp, "ERROR: Read beyond bounds of event: %d, run: %d\n", event->eventNo,
		    event->runNo);
	  }
	  length = -1; ERR = 1;
	  continue;
	}
	if(verbose){
	  if(logfp)
	    fprintf(logfp, "  Read %d channels of IU_MAM\n", reg_vals->nregs);
	}
	rawData++; length--;
	break;

      default:
	fprintf(stderr, "I don't recognize this component (0x%x)!\n",*rawData);

	if(logfp){
	  fprintf(logfp, "ERROR: Found a hardware component that is not known to RODD\n");
	  fprintf(logfp, "       Marker word: %d, Run: %d, Event: %d\n", 
		  *rawData, event->runNo, event->eventNo);
	}
	/* dumb test  */
	rawData++; length--;
	break;
      } /* end switch (*rawData) */
    
    }  /* end while length */
  } /* End of master loop over rawSize*/
  
  if (ERR == 1) {
    /* Something very bad took place */
   return(1);
  }
  
  /* at this point, we don't have any errors, write the data to the new event */
  
  if (lgd_adcs->nadc > 0) {
    size = sizeof_adc_values_t(lgd_adcs->nadc);
    data_removeGroup(event,GROUP_LGD_ADCS,0);
    group = data_addGroup(event, 100000, GROUP_LGD_ADCS, 0, size);
    memcpy(group, lgd_adcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of LGD_ADCS\n", lgd_adcs->nadc);
    }
  }
  /*  for now, we assume that the ONLY register is the MAM, so we have bypassed
       any need for dbmap lookup for the MAM  dsa */

  if (reg_vals->nregs > 0) {
     size = sizeof_reg_values_t(reg_vals->nregs);
     data_removeGroup(event,GROUP_MAM_REGS,0);
     group = data_addGroup(event, 100000, GROUP_MAM_REGS, 0, size);
     memcpy(group, reg_vals, size);
     if(verbose){
       if(logfp)
       fprintf(logfp, "   Wrote %d channels of MAM_REG\n", reg_vals->nregs);
     }
  }
  if (esum_adcs->nadc > 0) {
    size = sizeof_adc_values_t(esum_adcs->nadc);
    data_removeGroup(event,GROUP_ESUM_ADCS,0);
    group = data_addGroup(event, 100000, GROUP_ESUM_ADCS, 0, size);
    memcpy(group, esum_adcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of ESUM_ADCS\n", esum_adcs->nadc);
    }
  }
  if (rpd_adcs->nadc > 0) {
    size = sizeof_adc_values_t(rpd_adcs->nadc);
    data_removeGroup(event,GROUP_RPD_ADCS,0);
    group = data_addGroup(event, 100000, GROUP_RPD_ADCS, 0, size);
    memcpy(group, rpd_adcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of RPD_ADCS\n", rpd_adcs->nadc);
    }
  }
  if (cpv_adcs->nadc > 0) {
    size = sizeof_adc_values_t(cpv_adcs->nadc);
    data_removeGroup(event,GROUP_CPV_ADCS,0);
    group = data_addGroup(event, 100000, GROUP_CPV_ADCS, 0, size);
    memcpy(group, cpv_adcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of CPV_ADCS\n", cpv_adcs->nadc);
    }
  }
  if (upv_adcs->nadc > 0) {
    size = sizeof_adc_values_t(upv_adcs->nadc);
    data_removeGroup(event,GROUP_UPV_ADCS,0);
    group = data_addGroup(event, 100000, GROUP_UPV_ADCS, 0, size);
    memcpy(group, upv_adcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of UPV_ADCS\n", upv_adcs->nadc);
    }
  }
  if (bsd_adcs->nadc > 0) {
    size = sizeof_adc_values_t(bsd_adcs->nadc);
    data_removeGroup(event,GROUP_BSD_ADCS,0);
    group = data_addGroup(event, 100000, GROUP_BSD_ADCS, 0, size);
    memcpy(group, bsd_adcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of BSD_ADCS\n", bsd_adcs->nadc);
    }
  }
  if (bgv_adcs->nadc > 0) {
    size = sizeof_adc_values_t(bgv_adcs->nadc);
    data_removeGroup(event,GROUP_BGV_ADCS,0);
    group = data_addGroup(event, 100000, GROUP_BGV_ADCS, 0, size);
    memcpy(group, bgv_adcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of BGV_ADCS\n", bgv_adcs->nadc);
    }
  }
  if (tagger_tdcs->ntdc > 0) {
    size = sizeof_tdc_values_t(tagger_tdcs->ntdc);
    data_removeGroup(event,GROUP_TAGGER_TDCS,0);
    group = data_addGroup(event, 100000, GROUP_TAGGER_TDCS, 0, size);
    memcpy(group, tagger_tdcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of TAGGER_TDCS\n", tagger_tdcs->ntdc);
    }
  }
  if (rpd_tdcs->ntdc > 0) {
    size = sizeof_tdc_values_t(rpd_tdcs->ntdc);
    data_removeGroup(event,GROUP_RPD_TDCS,0);
    group = data_addGroup(event, 100000, GROUP_RPD_TDCS, 0, size);
    memcpy(group, rpd_tdcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of RPD_TDCS\n", rpd_tdcs->ntdc);
    }
  }
  if (cpv_tdcs->ntdc > 0) {
    size = sizeof_tdc_values_t(cpv_tdcs->ntdc);
    data_removeGroup(event,GROUP_CPV_TDCS,0);
    group = data_addGroup(event, 100000, GROUP_CPV_TDCS, 0, size);
    memcpy(group, cpv_tdcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of CPV_TDCS\n", cpv_tdcs->ntdc);
    }
  }
  if (upv_tdcs->ntdc > 0) {
    size = sizeof_tdc_values_t(upv_tdcs->ntdc);
    data_removeGroup(event,GROUP_UPV_TDCS,0);
    group = data_addGroup(event, 100000, GROUP_UPV_TDCS, 0, size);
    memcpy(group, upv_tdcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of UPV_TDCS\n", upv_tdcs->ntdc);
    }
  }
  if (bsd_tdcs_long->ntdc > 0) {
    size = sizeof_tdc_values_t(bsd_tdcs_long->ntdc);
    data_removeGroup(event,GROUP_BSD_TDCS_LONG,0);
    group = data_addGroup(event, 100000, GROUP_BSD_TDCS_LONG, 0, size);
    memcpy(group, bsd_tdcs_long, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of BSD_TDCS (long range)\n", bsd_tdcs_long->ntdc);
    }
  }
  if (bsd_tdcs_short->ntdc > 0) {
    size = sizeof_tdc_values_t(bsd_tdcs_short->ntdc);
    data_removeGroup(event,GROUP_BSD_TDCS_SHORT,0);
    group = data_addGroup(event, 100000, GROUP_BSD_TDCS_SHORT, 0, size);
    memcpy(group, bsd_tdcs_short, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of BSD_TDCS (short range)\n", bsd_tdcs_short->ntdc);
    }
  }
  if (bgv_tdcs->ntdc > 0) {
    size = sizeof_tdc_values_t(bgv_tdcs->ntdc);
    data_removeGroup(event,GROUP_BGV_TDCS,0);
    group = data_addGroup(event, 100000, GROUP_BGV_TDCS, 0, size);
    memcpy(group, bgv_tdcs, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of BGV_TDCS\n", bgv_tdcs->ntdc);
    }
  }
  if (scalers->nscalers > 0) {
    size = sizeof_scaler_values_t(scalers->nscalers);
    data_removeGroup(event,GROUP_SCALERS,0);
    group = data_addGroup(event, 100000, GROUP_SCALERS, 0, size);
    memcpy(group, scalers, size);
    if(verbose){
      if(logfp)
	fprintf(logfp, "   Wrote %d channels of SCALERS\n", scalers->nscalers);
    }
  }
  data_clean(event);  
  return(0);
}
static void reset_groups(void)
{
  esum_adcs->nadc=rpd_adcs->nadc = adc_vals->nadc = cpv_adcs->nadc = lgd_adcs->nadc = 0;
  tdc_vals->ntdc = rpd_tdcs->ntdc = cpv_tdcs->ntdc = tagger_tdcs->ntdc = 0;
  upv_adcs->nadc = upv_tdcs->ntdc = 0; 
  bsd_adcs->nadc = bsd_tdcs_short->ntdc = bsd_tdcs_long->ntdc = 0; 
  bgv_adcs->nadc = bgv_tdcs->ntdc = 0; 
  scaler_vals->nscalers = scalers->nscalers= 0;
  reg_vals->nregs = 0;
}
static void init_groups(void) 
{
  if(!adc_vals)
    if ((adc_vals = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL)
      malloc_error(__FILE__, __LINE__);      
  if(!tdc_vals)
    if ((tdc_vals = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__, __LINE__);
  if(!rpd_adcs)
    if ((rpd_adcs = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!lgd_adcs)
    if ((lgd_adcs = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL)
      malloc_error(__FILE__, __LINE__);
  if(!esum_adcs)
    if ((esum_adcs = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL)
      malloc_error(__FILE__, __LINE__);
  if(!cpv_adcs)
    if ((cpv_adcs = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL) 
      malloc_error(__FILE__, __LINE__);
  if(!upv_adcs)
    if ((upv_adcs = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL) 
      malloc_error(__FILE__, __LINE__);
  if(!bsd_adcs)
    if ((bsd_adcs = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL) 
      malloc_error(__FILE__, __LINE__);
  if(!bgv_adcs)
    if ((bgv_adcs = (adc_values_t *)malloc(sizeof_adc_values_t(1024))) == NULL) 
      malloc_error(__FILE__, __LINE__);
  if(!rpd_tdcs)
    if ((rpd_tdcs = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!cpv_tdcs)
    if ((cpv_tdcs = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!upv_tdcs)
    if ((upv_tdcs = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!bsd_tdcs_long)
    if ((bsd_tdcs_long = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!bsd_tdcs_short)
    if ((bsd_tdcs_short = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!bgv_tdcs)
    if ((bgv_tdcs = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!tagger_tdcs)
    if ((tagger_tdcs = (tdc_values_t *)malloc(sizeof_tdc_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!scaler_vals)
    if ((scaler_vals = (scaler_values_t *)malloc(sizeof_scaler_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!scalers)
    if ((scalers = (scaler_values_t *)malloc(sizeof_scaler_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  if(!reg_vals)
    if ((reg_vals = (reg_values_t *)malloc(sizeof_reg_values_t(1024))) == NULL)
      malloc_error(__FILE__,__LINE__);
  reset_groups();
}
static void malloc_error( char* source, int line )
{
    fprintf(stderr, "\n\nError in memory allocation!\n");
    fprintf(stderr, "Memory allocation failed at line %d in source file %s\n", line, source);
    fprintf(stderr, "Unable to continue, bailing out.\n");
    if(logfp){
      fprintf(logfp, "RODD terminated due to malloc error\n");
      fprintf(logfp, "----------------------------------------------\n");
      fclose(logfp);
    }
    exit(EXIT_FAILURE);
}
static void umap_error(int ROC, int slot, int channel) 
{
  fprintf(stderr, "Error in getting hardware mapping information, check log\n");
  fprintf(stderr, "   ROC: %d, Slot: %d, Channel %d\n", ROC, slot, channel);
  if(logfp){
    fprintf(logfp, "ERROR: The following did not map correctly: \n");
    fprintf(logfp, "   ROC: %d, Slot: %d, Channel %d\n", ROC, slot, channel);
    fprintf(logfp, "   Either incorrect slot to module mapping, or entry does not exist\n");
  }
}
static int dbmap(int ROCnum, adc_values_t* adc_group, 
		 tdc_values_t *tdc_group, scaler_values_t* scaler_group) 
{
  umapData_t umap;
  int i;

  if(scaler_group != NULL){
    for(i=0;i<scaler_group->nscalers;i++){
      if (umap_getInfo(ROCnum, scaler_group->scaler[i].slot, scaler_group->scaler[i].channel, &umap)) {
	umap_error(ROCnum, scaler_group->scaler[i].slot, scaler_group->scaler[i].channel);
	return 1;
      }
      if (umap.type != UMAP_SCALER) {
	umap_error(ROCnum, scaler_group->scaler[i].slot, scaler_group->scaler[i].channel);
	return 1;
      }
      if (umap.detector != UMAP_SCALER) continue; /* this channel is not used, ignore it */
      switch(umap.group) {
      case GROUP_SCALERS: 
	scalers->scaler[scalers->nscalers].channel = umap.swChannel;
	scalers->scaler[scalers->nscalers].value = scaler_group->scaler[i].value;
	scalers->nscalers++;
	break;
      default:
	/* we found a detector that hasn't been implemented, or an error */
	return 1;
      }
    }
  }

  if (adc_group != NULL) {
    for(i = 0; i < adc_group->nadc; i++) {
      if (umap_getInfo(ROCnum, adc_group->adc[i].slot, adc_group->adc[i].channel, &umap)) {
	umap_error(ROCnum, adc_group->adc[i].slot, adc_group->adc[i].channel);
	/*
	return 1;
	*/
      }
      if (umap.type != UMAP_ADC) {
	umap_error(ROCnum, adc_group->adc[i].slot, adc_group->adc[i].channel);
	/*
	return 1;
	*/
      }
      if (umap.detector == UMAP_UNUSED) continue; /* this channel is not used, ignore it */
      switch(umap.group) {
      case GROUP_RPD_ADCS: 
	rpd_adcs->adc[rpd_adcs->nadc].channel = umap.swChannel;
	rpd_adcs->adc[rpd_adcs->nadc].value = adc_group->adc[i].value;
	rpd_adcs->nadc++;
	break;
      case GROUP_LGD_ADCS: 
	lgd_adcs->adc[lgd_adcs->nadc].channel = umap.swChannel;
	lgd_adcs->adc[lgd_adcs->nadc].value = adc_group->adc[i].value;
	lgd_adcs->adc[lgd_adcs->nadc].mam = adc_group->adc[i].mam;
	lgd_adcs->adc[lgd_adcs->nadc].dor = adc_group->adc[i].dor;
	lgd_adcs->nadc++;
	break;
      case GROUP_CPV_ADCS: 
	cpv_adcs->adc[cpv_adcs->nadc].channel = umap.swChannel;
	cpv_adcs->adc[cpv_adcs->nadc].value = adc_group->adc[i].value;
	cpv_adcs->nadc++;
	break;
      case GROUP_UPV_ADCS: 
	upv_adcs->adc[upv_adcs->nadc].channel = umap.swChannel;
	upv_adcs->adc[upv_adcs->nadc].value = adc_group->adc[i].value;
	upv_adcs->nadc++;
	break;
      case GROUP_BSD_ADCS: 
	bsd_adcs->adc[bsd_adcs->nadc].channel = umap.swChannel;
	bsd_adcs->adc[bsd_adcs->nadc].value = adc_group->adc[i].value;
	bsd_adcs->nadc++;
	break;
      case GROUP_BGV_ADCS: 
	bgv_adcs->adc[bgv_adcs->nadc].channel = umap.swChannel;
	bgv_adcs->adc[bgv_adcs->nadc].value = adc_group->adc[i].value;
	bgv_adcs->nadc++;
	break;
      case GROUP_ESUM_ADCS:
	esum_adcs->adc[esum_adcs->nadc].channel = umap.swChannel;
	esum_adcs->adc[esum_adcs->nadc].value = adc_group->adc[i].value;
	esum_adcs->nadc++;
	break;
      default:
	/* we found a detector that hasn't been implemented, or an error */
	return 1;
      }
    }
  }
  
  if (tdc_group != NULL) {
    for(i = 0; i < tdc_group->ntdc; i++) {
      if (umap_getInfo(ROCnum, tdc_group->tdc[i].slot, tdc_group->tdc[i].channel, &umap)) {
	umap_error(ROCnum, tdc_group->tdc[i].slot, tdc_group->tdc[i].channel);
	return 1;
      }
      if (umap.type != UMAP_TDC) {
	umap_error(ROCnum, tdc_group->tdc[i].slot, tdc_group->tdc[i].channel);
	return 1;
      }
      if (umap.detector == UMAP_UNUSED) continue;   /* this channel is not used */
      switch(umap.group) {
      case GROUP_RPD_TDCS: 
	rpd_tdcs->tdc[rpd_tdcs->ntdc].channel = umap.swChannel;
	rpd_tdcs->tdc[rpd_tdcs->ntdc].le = tdc_group->tdc[i].le;
	rpd_tdcs->tdc[rpd_tdcs->ntdc].te = tdc_group->tdc[i].te;
	rpd_tdcs->ntdc++;
	break;
      case GROUP_TAGGER_TDCS: 
	tagger_tdcs->tdc[tagger_tdcs->ntdc].channel = umap.swChannel;
	tagger_tdcs->tdc[tagger_tdcs->ntdc].le = tdc_group->tdc[i].le;
	tagger_tdcs->tdc[tagger_tdcs->ntdc].te = tdc_group->tdc[i].te;
	tagger_tdcs->ntdc++;
	break;
      case GROUP_CPV_TDCS:
	cpv_tdcs->tdc[cpv_tdcs->ntdc].channel = umap.swChannel;
	cpv_tdcs->tdc[cpv_tdcs->ntdc].le = tdc_group->tdc[i].le;
	cpv_tdcs->tdc[cpv_tdcs->ntdc].te = tdc_group->tdc[i].te;
	cpv_tdcs->ntdc++;
	break;
      case GROUP_UPV_TDCS:
	upv_tdcs->tdc[upv_tdcs->ntdc].channel = umap.swChannel;
	upv_tdcs->tdc[upv_tdcs->ntdc].le = tdc_group->tdc[i].le;
	upv_tdcs->tdc[upv_tdcs->ntdc].te = tdc_group->tdc[i].te;
	upv_tdcs->ntdc++;
	break;
      case GROUP_BSD_TDCS_LONG:
	bsd_tdcs_long->tdc[bsd_tdcs_long->ntdc].channel = umap.swChannel;
	bsd_tdcs_long->tdc[bsd_tdcs_long->ntdc].le = tdc_group->tdc[i].le;
	bsd_tdcs_long->tdc[bsd_tdcs_long->ntdc].te = tdc_group->tdc[i].te;
	bsd_tdcs_long->ntdc++;
	break;
      case GROUP_BSD_TDCS_SHORT:
	bsd_tdcs_short->tdc[bsd_tdcs_short->ntdc].channel = umap.swChannel;
	bsd_tdcs_short->tdc[bsd_tdcs_short->ntdc].le = tdc_group->tdc[i].le;
	bsd_tdcs_short->tdc[bsd_tdcs_short->ntdc].te = tdc_group->tdc[i].te;
	bsd_tdcs_short->ntdc++;
	break;
      case GROUP_BGV_TDCS:
	bgv_tdcs->tdc[bgv_tdcs->ntdc].channel = umap.swChannel;
	bgv_tdcs->tdc[bgv_tdcs->ntdc].le = tdc_group->tdc[i].le;
	bgv_tdcs->tdc[bgv_tdcs->ntdc].te = tdc_group->tdc[i].te;
	bgv_tdcs->ntdc++;
	break;
      default: 
	/* found a detector that hasn't been implemented for this module, or an error */
	return 1;
      }
    }
  }
    
  return 0;
}
