/*
 * itypes.h - unpacked data - types (after the unpacker)
*/

#ifndef ITYPES_INCLUDED
#define ITYPES_INCLUDED


#include <ntypes.h>      /* for int32, etc... definitions */

typedef struct
{
  uint16 mam;
  uint16 dor;
  uint32 slot;
  uint32 channel;
  uint32 value;

} adc_value_t;

typedef struct
{
  int32 nadc;
  adc_value_t adc[1];
} adc_values_t;

#define sizeof_adc_values_t(N) \
       (sizeof(adc_values_t) + (((N)>0)?(N)-1:0)*sizeof(adc_value_t))

typedef struct
{
  uint32 slot;      
  uint32 channel;
  int32  le;         /* leading edge  */
  int32  te;         /* trailing edge */
} tdc_value_t;

typedef struct
{
  int32 ntdc;
  tdc_value_t tdc[1];
} tdc_values_t;

#define sizeof_tdc_values_t(N) \
       (sizeof(tdc_values_t) + (((N)>0)?(N)-1:0)*sizeof(tdc_value_t))

typedef struct
{
  uint32 slot;
  uint32 channel;
  uint32 value;
} scaler_value_t;

typedef struct
{
  int32  nscalers;
  scaler_value_t scaler[1];
} scaler_values_t;

#define sizeof_scaler_values_t(N) \
       (sizeof(scaler_values_t) + (((N)>0)?(N)-1:0)*sizeof(scaler_value_t))

typedef struct
{
  int32  nlatches;
  uint32 latch[1];
} latch_values_t;

#define sizeof_latch_values_t(N) \
       (sizeof(latch_values_t) + (((N)>0)?(N)-1:0)*sizeof(uint32))

typedef struct
{
  int32   channel;
  int32   adc;
  float32 fvalue;
} adc_hit_t;

typedef struct
{
  int32 nhits;
  adc_hit_t hit[1];
} adc_hits_t;

#define sizeof_adc_hits_t(N) \
       (sizeof(adc_hits_t) + (((N)>0)?(N)-1:0)*sizeof(adc_hit_t))

typedef struct
{
  int32 index;
  uint32 contents;
} reg_value_t;

typedef struct
{
  int32 nregs;
  reg_value_t reg[1];
} reg_values_t;

#define sizeof_reg_values_t(N) \
       (sizeof(reg_values_t) + (((N)>0)?(N)-1:0)*sizeof(reg_value_t))

typedef struct
{
  int32 nvalues;
  uint16 value[1];
} frozen_groups_t;

#define sizeof_frozen_groups_t(N) \
       (sizeof(frozen_groups_t) + (((N)>0)?(N)-1:0)*sizeof(uint16))

#endif
/* end file */
