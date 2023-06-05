/*
  $Log: pedestalsToMap.c,v $
  Revision 1.3  1999/07/29 10:40:54  radphi
  Changes made by gevmage@jlabs1

  Revision 1.2  1999/07/04 17:33:50  radphi
  Changes made by gevmage@jlabs1
  clarifying file name conventions

  Revision 1.1.1.1  1998/06/09 01:02:20  radphi
  Initial revision by radphi@urs2

 * Revision 1.2  1997/07/28  11:52:39  radphi
 * Changes made by lfcrob@dustbunny
 * Added other remaining detectors
 *
 * Revision 1.1.1.1  1997/05/24  00:46:54  radphi
 * Initial revision by lfcrob@jlabs2
 *
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pedestals.h>
#include <lgdGeom.h>

static const char rcsid[]="$Id: pedestalsToMap.c,v 1.3 1999/07/29 10:40:54 radphi Exp $";

void main(int argc, char *argv[])

{
  int iarg;
  char *fname;
  FILE *fp;
  int i;
  char tmp[1024];
  int detector = -1;
  int runNo = -1;
  int nChannels = -1;
  float *ped=NULL;
  int *thresh=NULL;
  int chan;
  float tmpPed;
  int tmpThresh;
  int nRead=0;

  for(iarg=1;iarg<argc;iarg++){
    fname = argv[iarg];
    fp = fopen(fname,"r");
    if(!fp){
      fprintf(stderr,"Failed to open %s\n",fname);      
      exit(0);
    }
    fscanf(fp,"%s %d %d",tmp,&runNo,&nChannels);
    if(!strcmp(tmp,"PED_LGD")){
       detector = PED_LGD;
    }
    else if(!strcmp(tmp,"PED_UPV")){
       detector = PED_UPV;
    }
    else if(!strcmp(tmp,"PED_BSD")){
       detector = PED_BSD;
    }
    else if(!strcmp(tmp,"PED_BGV")){
       detector = PED_BGV;
    }
    else if(!strcmp(tmp,"PED_CPV")){
       detector = PED_CPV;
    }
    else if(!strcmp(tmp,"PED_ADCS_TOP")){
       detector = PED_ADCS_TOP;
    }
    else if(!strcmp(tmp,"PED_ADCS_BOTTOM")){
       detector = PED_ADCS_BOTTOM;
    }

    else{
      fprintf(stderr,"Unknown detector : %s\n",tmp);
      exit(1);
    }
    /* Now read the pedestals */
    
    ped = malloc(nChannels*sizeof(float));
    memset(ped,0,nChannels*sizeof(float));
    thresh = malloc(nChannels*sizeof(int));
    memset(thresh,0,nChannels*sizeof(int));
    
    while(fscanf(fp,"%d %f %d",&chan,&tmpPed,&tmpThresh) == 3){
      ped[chan] = tmpPed;
      thresh[chan] = tmpThresh;
      nRead++;
    }
    fprintf(stderr,"Read %d pedestals from %s\n",nRead,fname);
    if(pedestalSave(detector,runNo,nChannels,ped,thresh)){
      fprintf(stderr,"Failed to save pedestals!\n");
      exit(0);
    }
  }
}



