/*Bug fix 6/99 by Tom Bogue.  last channel now used*/
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

/* CPVsave writes hits in CPV into the structures defined in 
 * the file iitypes.h                                
*/             
int CPVsave(itape_header_t *event,int nvdimcpv,int nhdimcpv,
            int nhmax,int *index,time_list_t *time_list) {
            
  cpv_hits_t *cpv_hits=NULL;
  int times[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float energy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int listdim=30;
  float timelist[30*100]; /* hopefully no more than 100 hits per 
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
  char chset[]="CPV ";
  char cvnames[5][5];
  char *chdet;
  int nvdim=nvdimcpv;
  int nhdim=nhdimcpv;
  int itrs=0;
  int numvs=0;
  int nhits;
  int len=4;
  
  strcpy(cvnames[0],"CPSA");
  strcpy(cvnames[1],"CPSB");
  strcpy(cvnames[2],"CPSC");
  strcpy(cvnames[3],"CPSD");
  strcpy(cvnames[4],"CPSE");

/* Malloc memory for the arrays */
  if ((numbv==NULL) && ((numbv=malloc(nvdimcpv*nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for numbv in CPVsave\n");
   return 0;
  } 
  if ((hits==NULL) && ((hits=malloc(nhdimcpv*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for hits in CPVsave\n");
   return 0;
  } 
  if ((sorthits==NULL) && ((sorthits=malloc(nhdimcpv*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for sorthits in CPVsave\n");
   return 0;
  } 
  if ((itra==NULL) && ((itra=malloc(nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for itra in CPVsave\n");
   return 0;
  } 

/* Loop over the detectors */
  for (key=0;key<5;key++) {
/* read the hits for the given key */
   chdet=cvnames[key];
   gfhits_(chset, chdet, &nvdim, &nhdim, &nhmax, &itrs, 
           &numvs, itra, numbv, hits, &nhits, 
           len, len);
         
/* First it is necessary to sort all the hits such that hit time
 * is increasing down the list
*/
   sort(hits,sorthits,nhdimcpv,nhits);
 
/* Now that hits are aligned we can sort them as described */
   timeAnalysisCPV (numbv,energy,times,timelist,listdim,
                 sorthits,nhits,key,nvdimcpv,nhdimcpv);   
  }  
/* now times contains number of DISTINGUISHABLE hits in a given channel
 * and timelist (array of times elements for each channel) 
 * contains the front edges of the hits
*/
  nhits=0;
  for (i=0;i<listdim;i++) { /* find the number of the hits in the CPV */
   if (energy[i]>0) nhits++;
  }    
/* Calculate the size of the structure*/
  size=sizeof_cpv_hits_t(nhits); 

/* Initialize the new data group*/
  cpv_hits=data_addGroup(event,BUFSIZE,GROUP_CPV_HITS,0,size);
  cpv_hits->nhits=nhits;
  
/* Fill the structure */
  k=0;
  for (i=0;i<listdim;i++) {
   if (energy[i]>0) {
    cpv_hits->hit[k].channel=i;
    cpv_hits->hit[k].energy=energy[i];
    cpv_hits->hit[k].times=times[i];
    cpv_hits->hit[k].tindex=*index;
    k++;
    for (l=0;l<times[i];l++) {
     time_list->element[*index].le=timelist[i+listdim*l]*1e9;
     (*index)++; /* This will shift the tindex */
    }  
   }
  }
  return 1;
 }   
