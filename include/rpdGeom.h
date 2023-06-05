/*
  $Log: rpdGeom.h,v $
  Revision 1.2  1997/06/02 18:38:30  radphi
  Initial revision by lfcrob@dustbunny
  Added some stuff...
  rpdGeomGetAngle(), rpdGeomGetName() rpdGeomGetMask(), rpdGetGeomParam()

 * Revision 1.1  1997/04/08  21:13:27  radphi
 * Initial revision
 *
  */

#include <encoderHandle.h>
#include <ntypes.h>

#ifndef RPDGEOM_H_INCLUDED
#define RPDGEOM_H_INCLUDED

int rpdGeomSetup(int runNo);
int rpdLocateXporter(encoderData_t *encoders);
vector3_t rpdGeomGetPosition(void);
float rpdGeomGetAngle(int channel);
const char *rpdGeomGetName(int channel);
int rpdGeomGetMask(int channel);
void rpdGetGeomParam(int *nChannels,int *nEChannels, int *nGChannels);
#endif
