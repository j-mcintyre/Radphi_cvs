/* written 6/99 by Tom Bogue*/
/* Standard libraries*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>


/* Misc Radphi */

#include <itypes.h>
#include <iitypes.h>
#include <disData.h>
#include <particleType.h>

/* Radphi I/O */

#include <dataIO.h>
#include <disIO.h>
#include <tapeIO.h>
#include <tapeData.h>

#include <mctypes.h>

#include "outputEvents.h"

/* BGVsave writes hits in BGV into the structures defined in iitypes.h */

int BGVsave(itape_header_t *event,
            int nvdimbgv,int nhdimbgv,int nhmaxbgv,
            int *index,time_list_t *time_list) {
 bgv_hits_t *bgv_hits=NULL;
 int nophot=0;
 float energy=0;
 int size=0;
 int i=0;
 int num=1;
 int end;

 int channelmap[] = {-1,0,1,2,3,4,5,6,7,8,9,10,11,12,
                     13,14,15,16,17,18,19,20,21,22,23};
/* channelmap[i] is the channel connected to the ith copy of the 
   volume 'BGVS' */
 
/* These are variables for reading hits */ 
  static int *numbv=NULL; 
  static float *hits=NULL;
  static int *itra=NULL;
  char chset[]="BGV ";
  char chdet[]="BGVS";
  int nvdim=nvdimbgv;
  int nhdim=nhdimbgv;
  int nhmax=nhmaxbgv;
  int itrs=0;
  int numvs[]={0,0,0,0,0,0,0};
  int nhts;
  int len=4;

  int hitn,hiteln;
  
/* Malloc memory for the arrays */
  if ((numbv==NULL) && ((numbv=malloc(nvdimbgv*nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for numbv in BGVsave\n");
   return 0;
  } 
  if ((hits==NULL) && ((hits=malloc(nhdimbgv*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for hits in BGVsave\n");
   return 0;
  } 
  if ((itra==NULL) && ((itra=malloc(nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for itra in BGVsave\n");
   return 0;
  } 
    
/* Read the hits */
  gfhits_(chset, chdet, &nvdim, &nhdim, &nhmax, &itrs, 
          numvs, itra, numbv, hits, &nhts, len, len);
         
/* Calculate the size of the structure*/
  size=sizeof_bgv_hits_t(nhts);

/* Initialize the new data group*/
 bgv_hits=data_addGroup(event,BUFSIZE,GROUP_BGV_HITS,0,size);
 bgv_hits->nhits=nhts;
 
/* Now simply loop over the number of hits and fill the structure*/
 for(i=0;i<nhts;i++) {
 
/* Now lets calculate the number of the channel*/
  bgv_hits->hit[i].channel=channelmap[numbv[i*nvdimbgv]];
 
/* Now fill the energy.*/
  bgv_hits->hit[i].energy[0]=hits[10+i*nhdimbgv];
/* hit[10+i*nhdimbgv] is Eupstream for the ith hit */
  bgv_hits->hit[i].energy[1]=hits[11+i*nhdimbgv];
/* hit[11+...] is Edownstream */
  end=0;
  time_list->element[*index].le = extract_bgv_t0_(&hits[i*nhdimbgv],&end)*1e9;
  bgv_hits->hit[i].tindex[0] = (*index)++;
  bgv_hits->hit[i].times[0] = 1;
/* tUpstream is extracted from LaPlace transform */
  end=1;
  time_list->element[*index].le = extract_bgv_t0_(&hits[i*nhdimbgv],&end)*1e9;
  bgv_hits->hit[i].tindex[1] = (*index)++;
  bgv_hits->hit[i].times[1] = 1;
/* tDownstream is extracted from LaPlace transform */
 }
 return 1;
}
