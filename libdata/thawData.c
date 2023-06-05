
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <disData.h>
#include <makeHits.h>
#include <unpackData.h>

adc_values_t *adc_vals=NULL;
adc_values_t *rpd_adcs=NULL;
adc_values_t *lgd_adcs=NULL;
adc_values_t *cpv_adcs=NULL;
adc_values_t *upv_adcs=NULL;
adc_values_t *bsd_adcs=NULL;
adc_values_t *bgv_adcs=NULL;
adc_values_t *esum_adcs=NULL;
tdc_values_t *tdc_vals=NULL;
tdc_values_t *rpd_tdcs=NULL;
tdc_values_t *cpv_tdcs=NULL;
tdc_values_t *upv_tdcs=NULL;
tdc_values_t *bsd_tdcs=NULL;
tdc_values_t *bgv_tdcs=NULL;
tdc_values_t *tagger_tdcs=NULL;
reg_values_t *reg_vals=NULL;

/* 
 * Frozen groups are the essential information from the ADC, TDC and MAM
 * groups compacted to take up a minimum of space in a disk file.  The
 * structure of frozen_groups_t is a sequence of blocks of 16-bit values,
 * each of which is prefixed by a group type and a number of group elements.
 */

int thawData(itape_header_t *event, int bufsize)
{
  int i, n;
  int size;
  frozen_groups_t *icebox;

  if ((icebox = data_getGroup(event,GROUP_FROZEN,0)) == NULL) {
    fprintf(stderr,"thawData: no frozen groups!\n");
    return(1);
  }
  for (i=0;i<icebox->nvalues;) {
    uint16 type=icebox->value[i++];
    uint16 count=icebox->value[i++];
    switch(type) {
      case GROUP_BSD_ADCS: 
        data_removeGroup(event,GROUP_BSD_ADCS,0);
        size = sizeof_adc_values_t(count);
        bsd_adcs = data_addGroup(event, 100000, GROUP_BSD_ADCS, 0, size);
        bsd_adcs->nadc = 0;
        for (n=0;n<count;n++) {
          bsd_adcs->adc[n].mam = (icebox->value[i] & 0x8000) >> 15;
          bsd_adcs->adc[n].dor = (icebox->value[i] & 0x4000) >> 14;
          bsd_adcs->adc[n].slot = 0;
          bsd_adcs->adc[n].channel = (icebox->value[i++] & 0xfff);
          bsd_adcs->adc[n].value = (icebox->value[i++] & 0xfff);
          bsd_adcs->nadc++;
        }
        break;
      case GROUP_BSD_TDCS: 
        data_removeGroup(event,GROUP_BSD_TDCS,0);
        size = sizeof_tdc_values_t(count);
        bsd_tdcs = data_addGroup(event, 100000, GROUP_BSD_TDCS, 0, size);
        bsd_tdcs->ntdc = 0;
        for (n=0;n<count;n++) {
          int time=(icebox->value[i] & 0x3ff);
          int channel=(icebox->value[i++] & 0xfe00) >> 10;
          bsd_tdcs->tdc[n].slot = 0;
          bsd_tdcs->tdc[n].channel = channel;
          bsd_tdcs->tdc[n].le = time;
          bsd_tdcs->tdc[n].te = 0;
          bsd_tdcs->ntdc++;
        }
        break;
      case GROUP_BGV_ADCS: 
        data_removeGroup(event,GROUP_BGV_ADCS,0);
        size = sizeof_adc_values_t(count);
        bgv_adcs = data_addGroup(event, 100000, GROUP_BGV_ADCS, 0, size);
        bgv_adcs->nadc = 0;
        for (n=0;n<count;n++) {
          bgv_adcs->adc[n].mam = (icebox->value[i] & 0x8000) >> 15;
          bgv_adcs->adc[n].dor = (icebox->value[i] & 0x4000) >> 14;
          bgv_adcs->adc[n].slot = 0;
          bgv_adcs->adc[n].channel = (icebox->value[i++] & 0xfff);
          bgv_adcs->adc[n].value = (icebox->value[i++] & 0xfff);
          bgv_adcs->nadc++;
        }
        break;
      case GROUP_BGV_TDCS: 
        data_removeGroup(event,GROUP_BGV_TDCS,0);
        size = sizeof_tdc_values_t(count);
        bgv_tdcs = data_addGroup(event, 100000, GROUP_BGV_TDCS, 0, size);
        bgv_tdcs->ntdc = 0;
        for (n=0;n<count;n++) {
          int time=(icebox->value[i] & 0x3ff);
          int channel=(icebox->value[i++] & 0xfe00) >> 10;
          bgv_tdcs->tdc[n].slot = 0;
          bgv_tdcs->tdc[n].channel = channel;
          bgv_tdcs->tdc[n].le = time;
          bgv_tdcs->tdc[n].te = 0;
          bgv_tdcs->ntdc++;
        }
        break;
      case GROUP_UPV_ADCS: 
        data_removeGroup(event,GROUP_UPV_ADCS,0);
        size = sizeof_adc_values_t(count);
        upv_adcs = data_addGroup(event, 100000, GROUP_UPV_ADCS, 0, size);
        upv_adcs->nadc = 0;
        for (n=0;n<count;n++) {
          upv_adcs->adc[n].mam = (icebox->value[i] & 0x8000) >> 15;
          upv_adcs->adc[n].dor = (icebox->value[i] & 0x4000) >> 14;
          upv_adcs->adc[n].slot = 0;
          upv_adcs->adc[n].channel = (icebox->value[i++] & 0xfff);
          upv_adcs->adc[n].value = (icebox->value[i++] & 0xfff);
          upv_adcs->nadc++;
        }
        break;
      case GROUP_UPV_TDCS: 
        data_removeGroup(event,GROUP_UPV_TDCS,0);
        size = sizeof_tdc_values_t(count);
        upv_tdcs = data_addGroup(event, 100000, GROUP_UPV_TDCS, 0, size);
        upv_tdcs->ntdc = 0;
        for (n=0;n<count;n++) {
          int time=(icebox->value[i] & 0x3ff);
          int channel=(icebox->value[i++] & 0xfe00) >> 10;
          upv_tdcs->tdc[n].slot = 0;
          upv_tdcs->tdc[n].channel = channel;
          upv_tdcs->tdc[n].le = time;
          upv_tdcs->tdc[n].te = 0;
          upv_tdcs->ntdc++;
        }
        break;
      case GROUP_CPV_ADCS: 
        data_removeGroup(event,GROUP_CPV_ADCS,0);
        size = sizeof_adc_values_t(count);
        cpv_adcs = data_addGroup(event, 100000, GROUP_CPV_ADCS, 0, size);
        cpv_adcs->nadc = 0;
        for (n=0;n<count;n++) {
          cpv_adcs->adc[n].mam = (icebox->value[i] & 0x8000) >> 15;
          cpv_adcs->adc[n].dor = (icebox->value[i] & 0x4000) >> 14;
          cpv_adcs->adc[n].slot = 0;
          cpv_adcs->adc[n].channel = (icebox->value[i++] & 0xfff);
          cpv_adcs->adc[n].value = (icebox->value[i++] & 0xfff);
          cpv_adcs->nadc++;
        }
        break;
      case GROUP_CPV_TDCS: 
        data_removeGroup(event,GROUP_CPV_TDCS,0);
        size = sizeof_tdc_values_t(count);
        cpv_tdcs = data_addGroup(event, 100000, GROUP_CPV_TDCS, 0, size);
        cpv_tdcs->ntdc = 0;
        for (n=0;n<count;n++) {
          int time=(icebox->value[i] & 0x3ff);
          int channel=(icebox->value[i++] & 0xfe00) >> 10;
          cpv_tdcs->tdc[n].slot = 0;
          cpv_tdcs->tdc[n].channel = channel;
          cpv_tdcs->tdc[n].le = time;
          cpv_tdcs->tdc[n].te = 0;
          cpv_tdcs->ntdc++;
        }
        break;
      case GROUP_LGD_ADCS: 
        data_removeGroup(event,GROUP_LGD_ADCS,0);
        size = sizeof_adc_values_t(count);
        lgd_adcs = data_addGroup(event, 100000, GROUP_LGD_ADCS, 0, size);
        lgd_adcs->nadc = 0;
        for (n=0;n<count;n++) {
          lgd_adcs->adc[n].mam = (icebox->value[i] & 0x8000) >> 15;
          lgd_adcs->adc[n].dor = (icebox->value[i] & 0x4000) >> 14;
          lgd_adcs->adc[n].slot = 0;
          lgd_adcs->adc[n].channel = (icebox->value[i++] & 0xfff);
          lgd_adcs->adc[n].value = (icebox->value[i++] & 0xfff);
          lgd_adcs->nadc++;
        }
        break;
      case GROUP_MAM_REGS: 
        data_removeGroup(event,GROUP_MAM_REGS,0);
        n = icebox->value[i++];
        size = sizeof_reg_values_t(n);
        reg_vals = data_addGroup(event, 100000, GROUP_MAM_REGS, 0, size);
        reg_vals->nregs = n;
        for (n=0;n<reg_vals->nregs;n++) {
          reg_vals->reg[n].index = n;
          reg_vals->reg[n].contents = 0;
        }
        for (n=0;n<count;n++) {
          int index=icebox->value[i++];
          reg_vals->reg[index].contents = icebox->value[i++];
          reg_vals->reg[index].contents += (icebox->value[i++] << 16);
        }
        break;
      case GROUP_TAGGER_TDCS: 
        data_removeGroup(event,GROUP_TAGGER_TDCS,0);
        size = sizeof_tdc_values_t(count);
        tagger_tdcs = data_addGroup(event, 100000, GROUP_TAGGER_TDCS, 0, size);
        tagger_tdcs->ntdc = 0;
        for (n=0;n<count;n++) {
          int time=(icebox->value[i] & 0x3ff);
          int channel=(icebox->value[i++] & 0xfe00) >> 10;
          tagger_tdcs->tdc[n].slot = 0;
          tagger_tdcs->tdc[n].channel = channel;
          tagger_tdcs->tdc[n].le = time;
          tagger_tdcs->tdc[n].te = 0;
          tagger_tdcs->ntdc++;
        }
        break;
      default:
        fprintf(stderr,"thawData: group %d not implemented!\n",type);
	break;
    }
  }
  return(0);
}

int freezeData(itape_header_t *event, int bufsize)
{
  int i, n=0;
  uint16 *size;
  int bytes;
  frozen_groups_t *icebox, *freezer;

  icebox = malloc(bufsize);
  if (icebox == NULL) {
    fprintf(stderr,"freezeData: malloc failed for frozen group!\n");
    return(1);
  }
  if (bsd_adcs = data_getGroup(event,GROUP_BSD_ADCS,0)) {
    icebox->value[n++] = GROUP_BSD_ADCS;
    size = &(icebox->value[n++]);
    *size = 0;
    for (i=0;i<bsd_adcs->nadc;i++) {
      icebox->value[n++] = (bsd_adcs->adc[i].mam>0? 0x8000 : 0x0)
                  + (bsd_adcs->adc[i].dor>0? 0x4000 : 0x0)
                  + (bsd_adcs->adc[i].channel & 0xfff);
      icebox->value[n++] = bsd_adcs->adc[i].value;
      (*size)++;
    }
    data_removeGroup(event,GROUP_BSD_ADCS,0);
  }
  if (bsd_tdcs = data_getGroup(event,GROUP_BSD_TDCS,0)) {
    icebox->value[n++] = GROUP_BSD_TDCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<bsd_tdcs->ntdc;i++) {
      int time=bsd_tdcs->tdc[i].le;
      if ((time > TMIN_BSD) && (time < TMAX_BSD)) {
        icebox->value[n++] = (bsd_tdcs->tdc[i].channel << 10)
                    + (bsd_tdcs->tdc[i].le & 0x3ff);
        (*size)++;
      }
    }
    data_removeGroup(event,GROUP_BSD_TDCS,0);
  }
  if (bgv_adcs = data_getGroup(event,GROUP_BGV_ADCS,0)) {
    icebox->value[n++] = GROUP_BGV_ADCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<bgv_adcs->nadc;i++) {
      icebox->value[n++] = (bgv_adcs->adc[i].mam>0? 0x8000 : 0x0)
                  + (bgv_adcs->adc[i].dor>0? 0x4000 : 0x0)
                  + (bgv_adcs->adc[i].channel & 0xfff);
      icebox->value[n++] = bgv_adcs->adc[i].value;
      (*size)++;
    }
    data_removeGroup(event,GROUP_BGV_ADCS,0);
  }
  if (bgv_tdcs = data_getGroup(event,GROUP_BGV_TDCS,0)) {
    icebox->value[n++] = GROUP_BGV_TDCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<bgv_tdcs->ntdc;i++) {
      int time=bgv_tdcs->tdc[i].le;
      if ((time > TMIN_BGV) && (time < TMAX_BGV)) {
        icebox->value[n++] = (bgv_tdcs->tdc[i].channel << 10)
                    + (bgv_tdcs->tdc[i].le & 0x3ff);
        (*size)++;
      }
    }
    data_removeGroup(event,GROUP_BGV_TDCS,0);
  }
  if (upv_adcs = data_getGroup(event,GROUP_UPV_ADCS,0)) {
    icebox->value[n++] = GROUP_UPV_ADCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<upv_adcs->nadc;i++) {
      icebox->value[n++] = (upv_adcs->adc[i].mam>0? 0x8000 : 0x0)
                  + (upv_adcs->adc[i].dor>0? 0x4000 : 0x0)
                  + (upv_adcs->adc[i].channel & 0xfff);
      icebox->value[n++] = upv_adcs->adc[i].value;
      (*size)++;
    }
    data_removeGroup(event,GROUP_UPV_ADCS,0);
  }
  if (upv_tdcs = data_getGroup(event,GROUP_UPV_TDCS,0)) {
    icebox->value[n++] = GROUP_UPV_TDCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<upv_tdcs->ntdc;i++) {
      int time=upv_tdcs->tdc[i].le;
      if ((time > TMIN_UPV) && (time < TMAX_UPV)) {
        icebox->value[n++] = (upv_tdcs->tdc[i].channel << 10)
                    + (upv_tdcs->tdc[i].le & 0x3ff);
        (*size)++;
      }
    }
    data_removeGroup(event,GROUP_UPV_TDCS,0);
  }
  if (cpv_adcs = data_getGroup(event,GROUP_CPV_ADCS,0)) {
    icebox->value[n++] = GROUP_CPV_ADCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<cpv_adcs->nadc;i++) {
      icebox->value[n++] = (cpv_adcs->adc[i].mam>0? 0x8000 : 0x0)
                  + (cpv_adcs->adc[i].dor>0? 0x4000 : 0x0)
                  + (cpv_adcs->adc[i].channel & 0xfff);
      icebox->value[n++] = cpv_adcs->adc[i].value;
      (*size)++;
    }
    data_removeGroup(event,GROUP_CPV_ADCS,0);
  }
  if (cpv_tdcs = data_getGroup(event,GROUP_CPV_TDCS,0)) {
    icebox->value[n++] = GROUP_CPV_TDCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<cpv_tdcs->ntdc;i++) {
      int time=cpv_tdcs->tdc[i].le;
      if ((time > TMIN_CPV) && (time < TMAX_CPV)) {
        icebox->value[n++] = (cpv_tdcs->tdc[i].channel << 10)
                    + (cpv_tdcs->tdc[i].le & 0x3ff);
        (*size)++;
      }
    }
    data_removeGroup(event,GROUP_CPV_TDCS,0);
  }
  if (lgd_adcs = data_getGroup(event,GROUP_LGD_ADCS,0)) {
    icebox->value[n++] = GROUP_LGD_ADCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<lgd_adcs->nadc;i++) {
      icebox->value[n++] = (lgd_adcs->adc[i].mam>0? 0x8000 : 0x0)
                  + (lgd_adcs->adc[i].dor>0? 0x4000 : 0x0)
                  + (lgd_adcs->adc[i].channel & 0xfff);
      icebox->value[n++] = lgd_adcs->adc[i].value;
      (*size)++;
    }
    data_removeGroup(event,GROUP_LGD_ADCS,0);
  }
  if (reg_vals = data_getGroup(event,GROUP_MAM_REGS,0)) {
    icebox->value[n++] = GROUP_MAM_REGS;
    size = &icebox->value[n++];
    icebox->value[n++] = reg_vals->nregs;
    *size = 0;
    for (i=0;i<reg_vals->nregs;i++) {
      if (reg_vals->reg[i].contents != 0) {
        icebox->value[n++] = reg_vals->reg[i].index;
        icebox->value[n++] = (reg_vals->reg[i].contents & 0xffff);
        icebox->value[n++] = (reg_vals->reg[i].contents >> 16);
        (*size)++;
      }
    }
    data_removeGroup(event,GROUP_MAM_REGS,0);
  }
  if (tagger_tdcs = data_getGroup(event,GROUP_TAGGER_TDCS,0)) {
    icebox->value[n++] = GROUP_TAGGER_TDCS;
    size = &icebox->value[n++];
    *size = 0;
    for (i=0;i<tagger_tdcs->ntdc;i++) {
      int time=tagger_tdcs->tdc[i].le;
      if ((time > TMIN_TAG) && (time < TMAX_TAG)) {
        icebox->value[n++] = (tagger_tdcs->tdc[i].channel << 10)
                    + (tagger_tdcs->tdc[i].le & 0x3ff);
        (*size)++;
      }
    }
    data_removeGroup(event,GROUP_TAGGER_TDCS,0);
  }
  icebox->nvalues = n;
  data_removeGroup(event,GROUP_FROZEN,0);
  bytes = sizeof_frozen_groups_t(n);
  freezer = data_addGroup(event, 100000, GROUP_FROZEN, 0, bytes);
  memcpy(freezer,icebox,bytes);
  free(icebox);
  return(0);
}
