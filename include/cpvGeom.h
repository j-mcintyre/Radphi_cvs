/*
  $Log: cpvGeom.h,v $
  Revision 1.1  1997/06/02 18:37:27  radphi
  Initial revision by lfcrob@dustbunny
  CPV geometry routines... bare bones

  */


#ifndef CVPGEOM_H_INCLUDED
#define CVPGEOM_H_INCLUDED

int cpvGeomSetup(int runNo);
void cpvGetGeomParam(int *channels,int *nHChannels, int *nVChannels);
const char *cpvGeomGetName(int channel);

#endif
