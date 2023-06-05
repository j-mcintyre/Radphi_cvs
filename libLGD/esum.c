/*
  $Log: esum.c,v $
  Revision 1.1  1997/06/02 18:39:37  radphi
  Initial revision by lfcrob@dustbunny
  Basic esum stuff

  */


#include <stdio.h>

#include <param.h>
#include <esum.h>

static const char rcsid[]="$Id: esum.c,v 1.1 1997/06/02 18:39:37 radphi Exp $";


static int nChannels;
static int totalChannel;



int esumSetup(int runNo)
{  
  if(param_loadDatabase("esum",runNo)){
    return(1);
  }
  if((nChannels = param_getInt("esum.nChannels",-1)) == -1)
    return(1);
  if((totalChannel = param_getInt("esum.Total",-1)) == -1)
    return(1);
  return(0);
}

void esumGetParam(int *channels,int *total)
{
  *channels = nChannels;
  *total = totalChannel;
}
