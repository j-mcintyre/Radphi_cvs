/* libdecode.c */

/*  
 *  These are the routines to decode the hardware modules, ADCs and TDCs.
 *  They decode one module at a time and increment the incoming data pointer.
 *  The struct definitions are not used, but are there for documentation.
 *  The geographical address or slot is returned by the functions.
 *  The model number is stored in the function name.
 */

/*   add decoding of discriminator bits for IU adc   dsa  */
/*   add scaler 1151 decoding   dsa   */
/*   add IUADC decoding   dsa/rtj  28/3/98 */
/*   add 1877 decoding   dsa  31/3/98 */
/*   add 1875 decoding   dsa  25/06/98 */


#include "itypes.h"
#include "ntypes.h"
#include <stdio.h>
#include "libdecode.h"

static int i, temp, chan, row, column;

typedef struct {
  unsigned slot    : 5;  /* MSBs */
  unsigned         : 4;
  unsigned range   : 1;
  unsigned channel : 7;
  unsigned         : 4;
  unsigned value   : 12; /* LSBs */
} data_word_1885_t;

typedef struct {
  unsigned ID_Reg      : 8;   /* MSBs */
  unsigned bit         : 1;
  unsigned slot        : 5;
  unsigned numChannels : 6;
  unsigned event       : 10;  /* LSBs */
} header_word_10C6_t;

typedef struct {
  unsigned zero      : 1;     /* MSBs */
  unsigned channel_a : 5;
  unsigned data_a    : 10;
  unsigned X         : 1;
  unsigned channel_b : 5;
  unsigned data_b    : 10;    /* LSBs */
} data_word_10C6_t;

typedef struct {
  unsigned channelID : 16;     /* MSBs */
  unsigned           : 1;
  unsigned MAMdiscri : 1;
  unsigned DORdiscri : 1;
  unsigned           : 1;
  unsigned value     : 12;    /* LSBs */
} data_word_IUadc_t;

typedef struct {
  unsigned slot      : 5;     /* MSBs */
  unsigned           : 11;
  unsigned parity    : 1;
  unsigned           : 1;
  unsigned buffnum   : 3;
  unsigned numhits   : 11;   /* LSBs */
} header_word_1877_t;

typedef struct {
  unsigned slot      : 5;     /* MSBs */
  unsigned parity    : 1;
  unsigned bnummod4  : 2;
  unsigned channel   : 7;
  unsigned sign      : 1;
  unsigned value     : 16;    /* LSBs */
} data_word_1877_t;

typedef struct {
  unsigned slot      : 5;     /* MSBs */
  unsigned event     : 3;
  unsigned range     : 1;
  unsigned           : 1;
  unsigned channel   : 6;
  unsigned           : 4;
  unsigned value     : 12;    /* LSBs */
} data_word_1875_t;


int decode_1885(int** Ev, int* length, adc_values_t *adc_Data)
{
  adc_Data->nadc = 96;                /* always reads out 96 channels */
  for (i = 0; i < adc_Data->nadc; i++) {
    temp = **Ev;
    adc_Data->adc[i].value = (temp & (0xfff));
    adc_Data->adc[i].channel = ((temp >> 16) & (0x7f));
    adc_Data->adc[i].slot = ((temp >> 27) & (0x1f));
    (*Ev)++; *length -= 1;
  }
  return 0;
}

int decode_1151(int** Ev, int* length, scaler_values_t *scaler_Data) 
{
  scaler_Data->nscalers = 16;                /* always reads out 16 channels */
  for (i = 0; i < scaler_Data->nscalers; i++) {
    temp = **Ev;
    scaler_Data->scaler[i].channel = i;
    scaler_Data->scaler[i].value = temp;
    scaler_Data->scaler[i].slot = 6;   /* be consistent with database  */
    (*Ev)++; *length -= 1;
  }
  return 0;
}

int decode_7200(int** Ev, int* length, scaler_values_t *scaler_Data) 
{
  scaler_Data->nscalers = 32;                /* always reads out 32 channels */
  for (i = 0; i < scaler_Data->nscalers; i++) {
    temp = **Ev;
    scaler_Data->scaler[i].value = temp;
    scaler_Data->scaler[i].channel = i; 
    scaler_Data->scaler[i].slot = 1; /* better be consistent with database */
    (*Ev)++; *length -= 1;
  }
  return 0;
}

int decode_block7200(int** Ev, int* length, scaler_values_t *scaler_Data) 
{
  scaler_Data->nscalers = 0;
  for (i = 0; i < 10000; i++) {
    temp = **Ev;
    if (temp == END_BLOCK_7200)
      break;
    scaler_Data->scaler[i].value = temp;
    scaler_Data->scaler[i].channel = i;
    scaler_Data->scaler[i].slot = 0;
    scaler_Data->nscalers++;
    (*Ev)++; *length -= 1;
  }
  return 0;
}

int decode_10C6(int** Ev, int* length, tdc_values_t *tdc_Data)
{

  uint32 header;

  typedef struct {   /* treat a data words as 2 16bit signed integers */
    int16 high;
    int16 low;
  } data_word_t;

  data_word_t *word;

  tdc_Data->ntdc = i = 0;
  while (**Ev != END_10C6) {

    header = **Ev;

  /*  
   *  stucture is one header word followed by data words
   *  header word will tell how many data words to expect next
   */

    if((header & 0x3ff) == 1) {
      tdc_Data->ntdc += ((header >> 10) & (0x3f));
      (*Ev)++; *length -= 1;
      for (; i < tdc_Data->ntdc; i++) {
        word = (void *)(*Ev);
        tdc_Data->tdc[i].slot = ((header >> 16) & (0x1f));
        tdc_Data->tdc[i].channel = ((word->high >> 10) & (0x1f));
        tdc_Data->tdc[i].le = (word->high & 0x3ff);
        /* bit 15 of low order may be 1 if that low word is a dummy */
        /* this happens when there is an odd number of channels */
        if (word->low > 0) {
          i++;
          tdc_Data->tdc[i].slot = ((header >> 16) & (0x1f));
  	  tdc_Data->tdc[i].channel = ((word->low >> 10) & (0x1f));
	  tdc_Data->tdc[i].le = (word->low & 0x3ff);
        }
        (*Ev)++; *length -= 1;
      }
    }
    else {
      (*Ev)++; *length -= 1; /* skip and data words before header */
    }
  }
  return 0;
}

int decode_IUmam(int** Ev, int* length, reg_values_t *mam_Data) 
{
  mam_Data->nregs = 10;                /* always reads out 10 registers */
  for (i = 0; i < mam_Data->nregs; i++) {
    temp = **Ev;
    mam_Data->reg[i].index = i;
    mam_Data->reg[i].contents = temp;
    (*Ev)++; *length -= 1;
  }
  return 0;
}

int decode_IUadc(int** Ev, int* length, adc_values_t *adc_Data) 
{
  temp = **Ev;
  adc_Data->nadc = 0;
  for (i=0; (temp != END_IU_ADC) && (*length > 0); i++){
    adc_Data->nadc++;
    adc_Data->adc[i].value = (temp & (0xfff));
    chan = (temp >> 16 & (0xffff));
    if (chan <= 0x1c1c){
      row = (chan >> 8 & 0xff);
      column = (chan & 0xff);
      adc_Data->adc[i].slot = 0; /* block read - assigned slot 0 */
      adc_Data->adc[i].channel = column + 28*row;  /* loading software channel */
    }
    else {
      adc_Data->adc[i].slot = (chan >> 8 & 0xff) - 0xE0; /* this really is a slot */
      adc_Data->adc[i].channel = (chan & 0xff);
    }        

    adc_Data->adc[i].dor = (temp >> 13 & 0x1);
    adc_Data->adc[i].mam = (temp >> 14 & 0x1);
    (*Ev)++; *length -= 1;
    temp = **Ev;
  }   /* end of for i loop */
  return 0;
}

int decode_1877(int** Ev, int* length, tdc_values_t *tdc_Data) 
{

  uint32 header;
  int time;


  tdc_Data->ntdc = 0;
  header = **Ev;

  if(header==FB_BUGGERED){
    /* fast bus time out on reading 1877 */
    return 1;
  }

  /* here we get the number of words and the slot from
     the  header; multiple hit tdc...  */

  (*Ev)++; *length -= 1;

  for (i=0; *length > 0; i++) {

    /* pull out the channel number from
       each data word
    */
    temp = **Ev;
    if (temp == END_1877){
      if(header==L2_FAIL_1877){
	/* Level 2 failed to reset 1877 in time */
	return(1);
      }
      break;
    }
    tdc_Data->tdc[tdc_Data->ntdc].slot = ((header >> 27) & (0x1f));
    tdc_Data->tdc[tdc_Data->ntdc].channel = (temp >> 17 & 0x7f) ;
    time = (temp & 0xffff);
    if ((i > 0) || (time != BOGUS_1877_CODE)) {
      /* In the case of a fast-clear inside the buffering interval
       * the 1877 puts out a junk first word that contains the value
       * in CSR 18 in its data portion and a channel number of
       * 0,2,4... (in decreasing frequency) in its channel portion.
       * The fast-clear is intended to work during this time, but
       * we discovered this extra word is a "feature" of the board.
       * Fortunately this data value is always equal to the lower
       * 16 bits of CSR 18 and always appears as the first word
       * after the header, so it is not hard to intercept.  The bad
       * news is that BOGUS_TDC_CODE needs to be changed if the value
       * written to CSR 18 is ever changed in the readout list. --RTJ 7/99
       */
      tdc_Data->tdc[tdc_Data->ntdc++].le = time;
    }
    (*Ev)++; *length -= 1;
  }
  return 0;
}

int decode_1875(int** Ev, int* length, tdc_values_t *tdc_Data) 
{
  tdc_Data->ntdc = 0;
  
  if((**Ev)==L2_FAIL_1875){
    /* Level 2 failed to reset 1875 in time */
    return 1;
  }
  if((**Ev)==FB_BUGGERED){
    /* fast bus time out on reading 1875 */
    return 1;
  }

  i=0; 
  while(((**Ev)!=END_1875)&&(*length>0)){
    /* pull out the channel number from
       each data word */
    temp = **Ev;
    tdc_Data->tdc[i].slot = ((temp & 0xf8000000) >> 27);
    tdc_Data->tdc[i].channel = ((temp & 0x003f0000) >> 16);
    tdc_Data->tdc[i].le = (temp & 0xfff);

    (*Ev)++;
    (tdc_Data->ntdc)++;
    *length -= 1;
    i++;
  }
  return 0;
}
/* 
 * $Log: libdecode.c,v $
 * Revision 1.14  2000/05/29 20:41:56  radphi
 * Changes made by radphi@urs3
 * unpacking of MAM and DOR (upper and lower disrims) for IU ADC
 * dsa
 *
 * Revision 1.13  1999/07/20 00:46:51  radphi
 * Changes made by jonesrt@jlabs1
 * added BOGUS_1877_CODE to libdecode.c -rtj
 *
 * Revision 1.12  1998/06/26 23:56:29  radphi
 * Changes made by radphi@jlabs1
 * added LRS 1875a decoding  D.S. Armstrong
 *
 * Revision 1.11  1998/04/01 20:29:46  radphi
 * Changes made by jonesrt@jlabs3
 * correct a bug for events with missing 10C6 header word - R.T.Jones
 *
 * Revision 1.10  1998/04/01  17:08:35  radphi
 * Changes made by radphi@jlabs2
 * fix bug in iuadc decoding and 10c6 decoding ds armstrong
 *
 * Revision 1.9  1998/03/31  21:23:42  radphi
 * Changes made by radphi@jlabs2
 * added slot to raw structures
 *
 * Revision 1.8  1998/03/31  20:00:46  radphi
 * Changes made by radphi@jlabs2
 * real decoding of 1877 TDC  D. Armstrong
 *
 * Revision 1.7  1998/03/31  04:14:40  radphi
 * Changes made by jonesrt@zeus
 * another typo - R.T.Jones
 *
 * Revision 1.6  1998/03/31  03:37:23  radphi
 * Changes made by jonesrt@zeus
 * getting the decoding working - R.T.Jones
 *
 * Revision 1.4  1998/03/30  21:36:19  radphi
 * Changes made by radphi@urs3
 * update iu adc unpacking, add upv
 *
 * Revision 1.3  1998/03/30  19:36:30  radphi
 * Changes made by jonesrt@
 * new implementation for March 1998 run period - R.T.Jones
 *
 * Revision 1.2  1997/05/20  09:07:01  radphi
 * Changes made by lfcrob@jlabs2
 * Added VME scalers
 *
 * Revision 1.1  1997/05/16  05:40:56  radphi
 * Initial revision by lfcrob@jlabs2
 * Module decoding routines - from rodd
 *
 * Revision 1.2  1997/05/14  20:06:09  radphi
 * Changes made by radphi@jlabs2
 * add scaler handling
 *
 * Revision 1.1.1.1  1997/04/22  14:36:42  radphi
 * Initial revision by lfcrob@dustbunny
 *
 * Revision 1.1.1.1  1997/04/22  14:26:49  radphi
 * Initial revision by lfcrob@dustbunny
 *
 * Revision 1.2  1997/03/22 18:23:32  toc
 * Keeps countdown of length
 * Modifies incoming data pointer as it walks down words
 *
 * Revision 1.1  1997/03/22 02:49:56  toc
 * Initial revision
 *
 */
