/* bug fix 6/99 by Tom Bogue.  Last channel should work now*/
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

/* RPDsave writes hitts in RPD into the structures defined in 
 * the file iitypes.h                                
*/             
int RPDsave(itape_header_t *event,int nvdimrpd,int nhdimrpd,
            int nhmax,int *index,time_list_t *time_list) {
  rpd_hits_t *rpd_hits=NULL;
  int times[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float energy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int listdim=18;
  float timelist[18*100]; /* hopefully no more than 100 hits per 
                            *channel per event
                            */ 
/* Actual numbers of channels in RPD */
  int chdim=12; 
  int channel[]={0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0,
                6, 7, 8, 9,10,11,12,13,14,15,16,17};     
        
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
  char chset[]="RPD ";
  char det1[]="DSEN";
  char det2[]="THSN";
  char *chdet;
  int nvdim=0;
  int nhdim=nhdimrpd;
  int itrs=0;
  int numvs=0;
  int nhits;
  int len=4;
 
/* Malloc memory for the arrays */
  if ((numbv==NULL) && ((numbv=malloc(nvdimrpd*nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for numbv in RPDsave\n");
   return 0;
  } 
  if ((hits==NULL) && ((hits=malloc(nhdimrpd*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for hits in RPDsave\n");
   return 0;
  } 
  if ((sorthits==NULL) && ((sorthits=malloc(nhdimrpd*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for sorthits in RPDsave\n");
   return 0;
  } 
  if ((itra==NULL) && ((itra=malloc(nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for itra in RPDsave\n");
   return 0;
  } 

/* Loop over the detectors */
  for (key=0;key<2;key++) {
/* read the hits for the given key */
   switch(key) {
    case 0 :
     chdet=det1;
     nvdim=nvdimrpd-1; /* for DSEN numbv if one-dimentional! */
     break;
    case 1 :
     chdet=det2;
     nvdim=nvdimrpd;
   }
   gfhits_(chset, chdet, &nvdim, &nhdim, &nhmax, &itrs, 
           &numvs, itra, numbv, hits, &nhits, 
           len, len);
         
/* First it is necessary to sort all the hits such that hit time
 * is increasing down the list
*/
   sort(hits,sorthits,nhdimrpd,nhits);
 
/* Now that hits are aligned we can sort them as described */
   timeAnalysis (channel,chdim,numbv,energy,times,timelist,listdim,
                 sorthits,nhits,key,nvdimrpd,nhdimrpd);   
  }  
/* now times contains number of DISTINGUISHABLE hits in a given channel
 * and timelist (array of times elements for each channel) 
 * contains the front edges of the hits
*/
  nhits=0;
  for (i=0;i<listdim;i++) { /* find the number of the hits in the RPD */
   if (energy[i]>0) nhits++;
  }    
/* Calculate the size of the structure*/
  size=sizeof_rpd_hits_t(nhits);

/* Initialize the new data group*/
  rpd_hits=data_addGroup(event,BUFSIZE,GROUP_RPD_HITS,0,size); 
  rpd_hits->nhits=nhits;
  
/* Fill the structure */
  k=0;
  for (i=0;i<listdim;i++) {
   if (energy[i]>0) {
    rpd_hits->hit[k].channel=i;
    rpd_hits->hit[k].energy=energy[i];
    rpd_hits->hit[k].times=times[i];
    rpd_hits->hit[k].tindex=*index;
    k++;
    for (l=0;l<times[i];l++) {
     time_list->element[*index].le=timelist[i+listdim*l]*1e9;
     (*index)++; /* This will shift the tindex */
    }  
   }
  }
  return 1;
 }   
