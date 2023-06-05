/* bug fix 6/99 by Tom Bogue to use last channel */
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
#include <lgdCluster.h>
#include <disData.h>
#include <particleType.h>

/* Radphi I/O */

#include <dataIO.h>
#include <disIO.h>
#include <tapeIO.h>
#include <tapeData.h>

#include <mctypes.h>

#include "outputEvents.h"

/* UPVsave writes hits in UPV into the structures defined in 
 * the file iitypes.h                                
*/             
int UPVsave(itape_header_t *event,int nvdimupv,int nhdimupv,
            int nhmax,int *index,time_list_t *time_list) {
  upv_hits_t *upv_hits=NULL;
  int times[]={0,0,0,0,0,0,0,0};
  float energy[]={0,0,0,0,0,0,0,0};
  int listdim=8;
  float timelist[8*100]; /* hopefully no more than 100 hits per 
                            *channel per event
                            */ 
/* the volume number is the channel number for the CPV */
        
  int size;
  int key;
  int i;
  int k;
  int l;
  int n;
 
 
/* These are variables for reading hits */ 
  static int *numbv=NULL; 
  static float *hits=NULL;
  static float *sorthits=NULL;
  static int *itra=NULL;
  char chset[]="UPV ";
  char chdet[]="UPVS";
  int nvdim=nvdimupv;
  int nhdim=nhdimupv;
  int itrs=0;
  int numvs=0;
  int nhits;
  int len=4;
 
/* Malloc memory for the arrays */
  if ((numbv==NULL) && ((numbv=malloc(nvdimupv*nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for numbv in UPVsave\n");
   return 0;
  } 
  if ((hits==NULL) && ((hits=malloc(nhdimupv*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for hits in UPVsave\n");
   return 0;
  } 
  if ((sorthits==NULL) && ((sorthits=malloc(nhdimupv*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for sorthits in UPVsave\n");
   return 0;
  } 
  if ((itra==NULL) && ((itra=malloc(nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for itra in UPVsave\n");
   return 0;
  } 
    
  key=0; /* In this casee just dummy index */
  gfhits_(chset, chdet, &nvdim, &nhdim, &nhmax, &itrs, 
          &numvs, itra, numbv, hits, &nhits, 
          len, len);
         
/* First it is necessary to sort all the hits such that hit time
 * is increasing down the list
*/
  sort(hits,sorthits,nhdimupv,nhits);
 
/* Now that hits are aligned we can sort them as described */
  timeAnalysisCPV (numbv,energy,times,timelist,listdim,
                   sorthits,nhits,key,nvdimupv,nhdimupv);   

/* now times contains number of DISTINGUISHABLE hits in a given channel
 * and timelist (array of times elements for each channel) 
 * contains the front edges of the hits
*/
  nhits=0;
  for (i=0;i<listdim;i++) { /* find the number of the hits in the UPV*/
   if (energy[i]>0) nhits++;
  }    
/* Calculate the size of the structure*/
  size=sizeof_upv_hits_t(nhits);

/* Initialize the new data group*/
  upv_hits=data_addGroup(event,BUFSIZE,GROUP_UPV_HITS,0,size); 
  upv_hits->nhits=nhits;
  
/* Fill the structure */
  k=0;
  for (i=0;i<listdim;i++) {
   if (energy[i]>0) {
    upv_hits->hit[k].channel=i;
    upv_hits->hit[k].energy=energy[i];
    upv_hits->hit[k].times=times[i];
    upv_hits->hit[k].tindex=*index;
    k++;
    for (l=0;l<times[i];l++) {
     time_list->element[*index].le=timelist[i+listdim*l]*1e9;
     (*index)++; /* This will shift the tindex */
    } 
   }
  } 
  return 1;
 }   
