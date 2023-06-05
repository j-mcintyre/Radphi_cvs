/* libdecode.h */

#include "itypes.h"

/* data stream marker words */
#define START_1885 0xda001885
#define END_1885 0xda005881

#define START_10C6 0xda0010c6
#define END_10C6 0xda006c01

#define START_1151 0xda001151
#define END_1151 0xda001511

#define START_NIM_7200 0xda0e7200
#define END_NIM_7200  0xda00027e

#define START_ECL_7200 0xda097200
#define END_ECL_7200  0xda000279

#define START_BLOCK_7200 0xda007200
#define END_BLOCK_7200  0xda000027

#define START_IU_ADC 0xda000adc
#define END_IU_ADC 0xda000cda

#define START_IU_MAM 0xda0001e2
#define END_IU_MAM 0xda002e10

#define START_1877 0xda001877
#define END_1877 0xda007781
#define L2_FAIL_1877 0xda18772f

#define START_1875 0xda001875
#define END_1875 0xda005781
#define L2_FAIL_1875 0xda18752f

#define FB_BUGGERED 0xfb000bad

#define BOGUS_1877_CODE 0x408

#define ADC_BOUNDARY 0xdaffffff



/*
 * these two are included for backwards compatability with
 * the first test data read out in Feb '97
 * when they are removed from here, they must be removed from parse.c
 */

#define START_1885_OLD 0xda0000aa 
#define END_1885_OLD 0xda0000bb

int decode_10C6(int**, int*, tdc_values_t*);
int decode_1885(int**, int*, adc_values_t*);
int decode_1151(int**, int*, scaler_values_t*);
int decode_7200(int**, int*, scaler_values_t*);
int decode_block7200(int**, int*, scaler_values_t*);
int decode_1877(int**, int*, tdc_values_t*);
int decode_1877_S(int**, int*, tdc_values_t*);
int decode_1875(int**, int*, tdc_values_t*);
int decode_IUadc(int**, int*, adc_values_t*);
int decode_IUmam(int**, int*, reg_values_t*);
