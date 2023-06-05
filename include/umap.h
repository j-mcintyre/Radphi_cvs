/*
  $Log: umap.h,v $
  Revision 1.7  1999/06/15 20:49:22  radphi

   updated to include BSD and BGV detectors as part of hardware map
      D.S. Armstrong

  Revision 1.6  1998/06/29 20:04:44  radphi
  Changes made by jonesrt@urs3
  updates for new level i and level ii unpacking - RTJ

 * Revision 1.5  1998/03/31  02:35:14  radphi
 * Changes made by jonesrt@zeus
 * new types for new readout modules - R.T.Jones
 *
 * Revision 1.4  1997/05/20  09:02:02  radphi
 * Changes made by lfcrob@jlabs2
 * ddd
 *
 * Revision 1.3  1997/05/20  04:34:48  radphi
 * Changes made by radphi@jlabs2
 * Added ESUM
 *
 * Revision 1.2  1997/05/14  05:44:18  radphi
 * Changes made by lfcrob@jlabs2
 * Added routine to map swchannel->hw channel
 *
 * Revision 1.1  1997/02/26  19:58:04  lfcrob
 * Initial revision
 *
  */

#ifndef UNPACK_H_INCLUDED
#define UNPACK_H_INCLUDED


#define UMAP_UNKNOWN 0x1

#define UMAP_UNUSED 0x2
#define UMAP_LGD 0x4
#define UMAP_SCALER 0x8
#define UMAP_RPD 0x10
#define UMAP_CPV 0x20
#define UMAP_TAGGER 0x40
#define UMAP_ESUM 0x80
#define UMAP_UPV 0x100
#define UMAP_BSD 0x200
#define UMAP_BGV 0x400

#define UMAP_VME 0x2
#define UMAP_FASTBUS 0x4
#define UMAP_CAMAC 0x8

#define UMAP_ADC 0x2
#define UMAP_TDC 0x4
#define UMAP_SCALER 0x8
#define UMAP_REG 0x10

#define UMAP_LRS_1885M 0x2
#define UMAP_STK_10C6 0x4
#define UMAP_LRS_1151 0x8
#define UMAP_STK_7200 0x10
#define UMAP_IU_ADC 0x20
#define UMAP_IU_MAM 0x40
#define UMAP_LRS_1877 0x80
#define UMAP_LRS_1875 0x100

typedef struct{
  int group;
  int detector;
  int swChannel;
} umapChannel_t;

typedef struct{
  int model;
  int moduleType;
  int bus;
  int slot;
  int nChannels;
  umapChannel_t *channel;
} umapModule_t;

typedef struct{
  int rocNum;
  int nModules;
  umapModule_t **module;  
} umapRoc_t;

typedef struct{
  int nRocs;
  umapRoc_t **roc;
} umap_t;

typedef struct{
  int swChannel;
  int detector;
  int group;
  int model;
  int type;
  int bus;
} umapData_t;

void umapDataPrint(FILE *fp,umapData_t *data);
int umap_loadDb(int runNo);
int umap_getInfo(int roc, int slot, int channel, umapData_t *data);
int umap_getHW(int detector, int module, int channel,int *slot,int *hwChannel);
#endif
