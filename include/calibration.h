/*
  $Log: calibration.h,v $
  Revision 1.4  2001/06/03 17:54:33  radphi
  Changes made by crsteffe@mantrid00
  has provisions for lgdtune constants

  Revision 1.3  2000/07/27 03:36:30  radphi
  Changes made by gevmage@urs3
  added proper updating of constants from map

  Revision 1.2  2000/06/18 06:37:07  radphi
  Changes made by jonesrt@grendl
  updates to ccLoad for bsd,bgv,upv groups -rtj

 * Revision 1.1  1997/06/02  18:36:59  radphi
 * Initial revision by lfcrob@dustbunny
 * Calibration constant routine headers
 *
  */

#ifndef CALIBRATION_H_INCLUDED
#define CALIBRATION_H_INCLUDED

#define CC_LGD 0x1
#define CC_RPD 0x2
#define CC_CPV 0x4
#define CC_ESUM 0x8
#define CC_UPV 0x10
#define CC_BGV 0x20
#define CC_BSD 0x40

#define PARAMETER_EPSILON 0x1
#define PARAMETER_BETA 0x2

int ccLoad(int detector, int runNo, int nChannels, float *cc, int *actualRun);
int ccSave(int detector, int runNo, int nChannels, float *cc);
int ccParameterLoad(int parameter_name, int runNo, 
		    float *parameter, int *actualRun);
int lgdTuneLoad(int runNo, int nChannels, float *cc, int *actualRun);

#endif 
