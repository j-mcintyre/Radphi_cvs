/*
  $Log: cpvGeom.c,v $
  Revision 1.1.1.1  1998/06/29 21:43:58  radphi
  Initial revision by jonesrt@urs3

 * Revision 1.1  1997/06/02  18:39:26  radphi
 * Initial revision by lfcrob@dustbunny
 * Basic cpv geometry stuff
 *
  */


#include <stdio.h>

#include <param.h>
#include <cpvGeom.h>

static const char rcsid[]="$Id: cpvGeom.c,v 1.1.1.1 1998/06/29 21:43:58 radphi Exp $";

static int nChannels;
static int H_channelZero;
static int V_channelZero;

int cpvGeomSetup(int runNo)
{
  vector3_t dft3v={-1,-1,-1};
  
  if(param_loadDatabase("cpvGeom",runNo)){
    return(1);
  }
  if((nChannels = param_getInt("cpvGeom.nChannels",-1)) == -1)
    return(1);
  if((H_channelZero = param_getInt("cpvGeom.H.channelZero",-1)) == -1)
    return(1);
  if((V_channelZero = param_getInt("cpvGeom.V.channelZero",-1)) == -1)
    return(1);
  return(0);
}

void cpvGetGeomParam(int *channels,int *hchannels, int *vchannels)
{
  *channels = nChannels;
  if(H_channelZero < V_channelZero){
    *hchannels = V_channelZero;
    *vchannels = nChannels-V_channelZero;
  }
  else{
    *hchannels = H_channelZero;
    *vchannels = nChannels-H_channelZero;
  }
}
const char *cpvGeomGetName(int channel)
{
  char tmp[80];

  sprintf(tmp,"cpvGeom.channel.%d.name",channel);
  return(param_getValue(tmp));
}
