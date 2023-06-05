/*
 * February 13, 2003
 *   Added ability to store the actual photoelectron count from explicit
 *   Cerenkov following instead of the model yield smeared by Poisson
 *   statistics.  To enable this behavior, compile with the following.
 * #define FOLLOW_CERENKOV true
 * - Richard Jones
 */

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
#include "digitizer.h"

/* LGDsave writes hits in LGD into the structures defined in 
 * the file lgdCluster.h
 *                                
 * This is the map of lgd wall
 *
 *			wal2 (1)
 *
 *		wal1(2)		wal1(1)
 *
 *			wal2(2)
 *
 * looking upstream
*/             
int LGDsave(itape_header_t *event,
            int nvdimlgd,int nhdimlgd,int nhmaxlgd) {
  lgd_hits_t *lgd_hits=NULL;
  int size;
  int i;
 
/* These are variables for reading hits */ 
  static int *numbv=NULL; 
  static float *hits=NULL;
  static int *itra=NULL;
  char chset[]="CALO";
  char chdet[]="BLOK"; 
  int nvdim=nvdimlgd;
  int nhdim=nhdimlgd;
  int nhmax=nhmaxlgd;
  int itrs=0;
  int numvs[]={0,0,0,0,0,0,0};
  int nhts;
  int len=4;
  const int nlgd=LGD_ROWS*LGD_COLUMNS;
  int *lgdsum;
  int channel;
  
/* Malloc memory for the arrays */
  if ((numbv==NULL) && ((numbv=malloc(nvdimlgd*nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for numbv in LGDsave\n");
   return 0;
  } 
  if ((hits==NULL) && ((hits=malloc(nhdimlgd*nhmax*sizeof(float)))==NULL)) {
   printf("Fail to accocate memory for hits in LGDsave\n");
   return 0;
  } 
  if ((itra==NULL) && ((itra=malloc(nhmax*sizeof(int)))==NULL)) {
   printf("Fail to accocate memory for itra in LGDsave\n");
   return 0;
  } 
  if ((lgdsum=calloc(nlgd,sizeof(float)))==NULL) {
   printf("Fail to accocate memory for sumlgd in LGDsave\n");
   return 0;
  } 

/* Read the hits */

#if defined FOLLOW_CERENKOV
  strcpy(chdet,"PMTC");
#endif
  gfhits_(chset, chdet, &nvdim, &nhdim, &nhmax, &itrs, 
          numvs, itra, numbv, hits, &nhts, 
          len, len);
         
/* Calculate the size of the structure*/
  size=sizeof_lgd_hits_t(nhts);

/* Initialize the new data group*/
  lgd_hits=data_addGroup(event,BUFSIZE,GROUP_LGD_HITS,0,size);
  lgd_hits->nhits=0;
 
/* Now simply loop over the number of hits and fill the structure */
  for(i=0;i<nhts;i++) {
   const int one=1;
   int block[2];
   int nophot;
#if defined FOLLOW_CERENKOV
  /* Here all we have to do is to decide what block was being watched
   * by this phototube, and assign the channel number.  The total
   * photoelectron yield is stored in the hits array.
   */
   block[0]=(hits[0+nhdimlgd*i]/4.0)+(LGD_ROWS/2.0)+1;
   block[1]=(hits[1+nhdimlgd*i]/4.0)+(LGD_COLUMNS/2.0)+1;
   nophot=hits[5+nhdimlgd*i];
#else
  /* Here we have to do two things: first, calculate the channel number
   * that corresponds to the given block based on its indices and
   * second, convert the readings of the hits array into units of GeV
   * after randomizing the initial value according to the Poisson
   * distribution.
   */
   if (numbv[0+nvdimlgd*i]>0) {
    block[1]=numbv[2+nvdimlgd*i];
    if (numbv[0+nvdimlgd*i]==2) {
     block[0]=numbv[1+nvdimlgd*i]; /* it is in wal1(2) */
    } 
    else {
     block[0]=numbv[1+nvdimlgd*i]+2+(LGD_ROWS-2)/2; /* it is in wal1(1) */
    } 
   }
   else {
    block[0]=numbv[4+nvdimlgd*i]+(LGD_ROWS-2)/2; /* it is in wal2 */
    if (numbv[3+nvdimlgd*i]==2) {
     block[1]=numbv[5+nvdimlgd*i];
    } 
    else {
     block[1]=numbv[5+nvdimlgd*i]+2+(LGD_COLUMNS-2)/2;
    }   
   }
  /* randomize photoelectron yield according to Poisson statistics */
   gpoiss_(&hits[5+nhdimlgd*i],&nophot,(int*)&one);
#endif
   channel=(block[1]-1)*LGD_ROWS+(block[0]-1);
   if (channel < 0 || channel >= nlgd) {
    printf("Warning in LGDsave: attempt to save a hit in LGD "
           "row %d, column %d, channel %d outside the physical "
           "limits of the detector! Hit will be dropped, "
           "moving on...", block[1], block[0], channel);
   }
   else {
    lgdsum[channel]+=nophot;
   }
  }

/* Now save the total photoelectron yield */
  for (channel=0; channel < nlgd; channel++) {
   int nophot=lgdsum[channel];
   if (nophot > 0) {
    int ihit=lgd_hits->nhits++;
    lgd_hits->hit[ihit].channel=channel;
    lgd_hits->hit[ihit].energy=nophot*LGD_GEV_PER_PE;
   }
  }
  free(lgdsum);
  return 1;
}
