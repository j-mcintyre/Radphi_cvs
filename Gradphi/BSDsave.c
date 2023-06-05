/* 6/99 by Tom Bogue. Bulk of this routine was copied from RPDhits*/
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


/* BSDsave writes hits in BSD into the structures defined in 
 * the file iitypes.h                                
*/             
int BSDsave(itape_header_t *event,int nvdimbsd,int nhdimbsd,
            int nhmax,int *index,time_list_t *time_list)
{
  bsd_hits_t *bsd_hits=NULL;
  int times[]={0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float energy[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int listdim=48;
  float timelist[48*100]; /* hopefully no more than 100 hits per 
                            *channel per event
                            */ 
/* Actual numbers of channels in BSD */
  int chdim=24; 
  int channel[]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,0,0,0,0,0,0,0,0,0,0,0,0,
		 12,13,14,15,16,17,18,19,20,21,22,23,0,0,0,0,0,0,0,0,0,0,0,0,
		 24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};
        
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
  char chset[]="BSD ";
  char det1[]="BSR1";
  char det2[]="BSR2";
  char det3[]="BSR3";
  char *chdet;
  int nvdim=nvdimbsd;
  int nhdim=nhdimbsd;
  int itrs=0;
  int numvs=0;
  int nhits;
  int len=4;

/* Malloc memory for the arrays */
  if ((numbv==NULL) && ((numbv=malloc(nvdimbsd*nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for numbv in BSDsave\n");
   return 0;
  } 
  if ((hits==NULL) && ((hits=malloc(nhdimbsd*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for hits in BSDsave\n");
   return 0;
  } 
  if ((sorthits==NULL) && ((sorthits=malloc(nhdimbsd*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for sorthits in BSDsave\n");
   return 0;
  } 
  if ((itra==NULL) && ((itra=malloc(nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for itra in BSDsave\n");
   return 0;
  } 
/* Loop over the detectors */
  for (key=0;key<3;key++) {
/* read the hits for the given key */
   switch(key) {
    case 0 :
     chdet=det1;
     break;
    case 1 :
     chdet=det2;
     break;
    case 2 :
     chdet=det3;
   }
   gfhits_(chset, chdet, &nvdim, &nhdim, &nhmax, &itrs, 
           &numvs, itra, numbv, hits, &nhits, 
           len, len);

/* First it is necessary to sort all the hits such that hit time
 * is increasing down the list
 */
   sort(hits,sorthits,nhdimbsd,nhits);

/* Now that hits are aligned we can sort them as described
 * since the individual scintillators in the BSD haven't been 
 * split yet, we need a special routine for the BSD
 */
   timeAnalysisBSD (channel,chdim,numbv,energy,times,timelist,listdim,
                    sorthits,nhits,key,nvdimbsd,nhdimbsd);
  }
/* now times contains number of DISTINGUISHABLE hits in a given channel
 * and timelist (array of times elements for each channel) 
 * contains the front edges of the hits
 */
  nhits=0;
  for (i=0;i<listdim;i++) { /* find the number of the hits in the BSD */
   if (energy[i]>0) nhits++;
  }
/* Calculate the size of the structure */
  size=sizeof_bsd_hits_t(nhits);

/* Initialize the new data group */
  bsd_hits=data_addGroup(event,BUFSIZE,GROUP_BSD_HITS,0,size);
  bsd_hits->nhits=nhits;

/* Fill the structure */
  k=0;
  for (i=0;i<listdim;i++) {
   if (energy[i]>0) {
    bsd_hits->hit[k].channel=i;
    bsd_hits->hit[k].energy=energy[i];
    bsd_hits->hit[k].times=times[i];
    bsd_hits->hit[k].tindex=*index;
    k++;
    for (l=0;l<times[i];l++) {
     time_list->element[*index].le=timelist[i+listdim*l]*1e9;
     (*index)++; /* This will shift the tindex */
    }  
   }
  }
  return 1;
}
