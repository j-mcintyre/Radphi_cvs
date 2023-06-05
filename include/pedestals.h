/*
  $Log: pedestals.h,v $
  Revision 1.5  2000/01/18 03:58:56  radphi
  Changes made by jonesrt@golem
  (1) Fixed endian-conversion package (libdata) for MonteCarlo itape files
  (2) Major upgrade to map manager package to support little-endian machines
  (3) Removed redundant reference (rpd/upv) in pedestals.h
  -rtj

  Revision 1.4  1999/07/09 02:55:33  radphi
  Changes made by radphi@urs2
  Removed redundant PED_CPV macro
  P. Rubin

  Revision 1.3  1999/07/08 20:50:27  radphi
  Changes made by radphi@urs1
  Returned pedestal values from old configuration to make backwards compatible
  P. Rubin

  Revision 1.2  1999/07/04 17:46:09  radphi
  Changes made by radphi@urs1
  has 1999 detectors, and adc control stuff

  Revision 1.1  1997/05/16 19:20:07  radphi
  Initial revision by lfcrob@jlabs2
  pedestal loading routines

  */

#ifndef PEDESTALS_H_INCLUDED
#define PEDESTALS_H_INCLUDED

#define PED_LGD 0x1
#define PED_UPV 0x2
#define PED_BSD 0x4
#define PED_BGV 0x8
#define PED_CPV 0x10
#define PED_ADCS_TOP 0x1000
#define PED_ADCS_BOTTOM 0x2000
#define PED_RPD 0x20
#define PED_ESUM 0x80

int pedestalLoad(int detector, int runNo, int nChannels, float *pedestal, int *threshold, int *actualRun);
int pedestalSave(int detector, int runNo, int nChannels, float *pedestal, int *threshold);
#endif 
