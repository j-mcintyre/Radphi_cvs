/*
* outputEvents: saves the event tracked by GEANT in the output file 
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
#include <eventType.h>

/* Radphi I/O */

#include <dataIO.h>
#include <disIO.h>
#include <tapeIO.h>
#include <tapeData.h>
#include <triggerType.h>

#include <mctypes.h>

#include "outputEvents.h"

int save_next_event_(int *runNo)

/* saves the hit in the corresponding structure (cpv, upv etc.)
 * it reads the hits, determines which detector was hit
 * and then just calls the corresponding subroutine to save it.
 * maxevent is the maximum number of events to be triggered in 
 * the current 'run'
*/ 
   
{
  static FILE *fpOut=NULL;
  static itape_header_t *event=NULL;
  static time_list_t *time_list=NULL; /* This buffer has an arbitrary large size, 
                                * to copied to real_time_list in the end
                               */ 
  static time_list_t *real_time_list=NULL; /* This one will have the actual size */
  
  mc_event_t *mc_event=NULL; /* real array */                          
  mc_vertex_t *mc_vertex=NULL;
  
/* This is the name of the output file*/   
  char filename[120]="simData.itape";
  
/* some other variables*/
  static int eventNo=0;
  static int spillNo=1;
  int ret;
  int i;
  int size;
  int nubuf=0;
  float ubuf[BUFSIZE];
  float vert[3];
  float pvert[3];
  int tindex=0;
    
/* These variables are to be used when fetching hits; see 
 * the description of gfhits in GEANT manual for details
*/
  int nvdimlgd=7;
  int nhdimlgd=6;
  int nvdimcpv=1;
  int nhdimcpv=7;
  int nvdimrpd=2;
  int nhdimrpd=7;
  int nvdimupv=1;
  int nhdimupv=7;
  int nvdimbsd=1;
  int nhdimbsd=7;
  int nvdimbgv=1;
  int nhdimbgv=12;
  int nhmax=800;
  int itrs;
  int len;
  int nhits;
  int numvs=0;
  
/* some dummy variables */
  int itra=1;
  int nvtx=1;
  int ipart;
  int ntbeam=0;
  int nttarg=0;
  float tofg=0;
  
/* Open the file for output; warn if error
 * the file is opened just once
 */
  if ((fpOut == NULL) && ((fpOut=fopen(filename,"w")) == NULL)) {
    fprintf(stderr,"WARNING - failed to open output event file %s\n",filename);
    return 0;
  }

/* Malloc the event buffer; do once */
  if ((event==NULL) && ((event=malloc(BUFSIZE)) == NULL)) {
    fprintf(stderr,"ERROR - failed to allocating output event buffer!\n");
    fclose(fpOut);
    return 0;
  }
  
  ++eventNo;
  
/* Initialize the event*/  
  data_newItape(event);
  event->eventType=EV_DATA;
  event->runNo = *runNo;
  event->spillNo = eventNo/100+1;
  event->eventNo = eventNo;
  event->trigger = TRIG_MC;
  event->latch = GRADPHI_VERSION_RELEASE;
  
/* Some arbitrary large number of hits for time_list */
  nhits=1000; /* Just some big number */  
 
/* Calculate the size of the time_event buffer and malloc it*/
  size=sizeof_time_list_t(nhits);
  if ((time_list==NULL) && ((time_list=malloc(size))==NULL)) {
   fprintf(stderr,"ERROR-fail to allocate memory for output time list!\n");
   fclose(fpOut);
   return 0;
  }
  
/* now fetch the user words array attached to the first event */
  {
   int nwrung=0;
   float irung[50];
   gfrung_(&nwrung,irung,&nubuf,ubuf);
  } 

/* recalculate size and add new data group to the event */
  size=nubuf*4;
  if (size>0) {
   mc_event=data_addGroup(event,BUFSIZE,GROUP_MC_EVENT,0,size);
   memcpy(mc_event,ubuf,size);
  } 
  
/* now fetch the user words array attached to the first vertex */
  gfvert_(&nvtx,vert,&ntbeam,&nttarg,&tofg,ubuf,&nubuf);
  
/* recalculate size and add new data group to the event */
  size=nubuf*4;
  if (size>0) {
   mc_vertex=data_addGroup(event,BUFSIZE,GROUP_MC_EVENT,1,size);
   memcpy(mc_vertex,ubuf,size);
  } 
  
/* Lets put a dummy size in the time_list (just for the time being, correct 
 * when adding the actual group
*/
  time_list->nelements=1;    

/* Now fetch the hits and call the corresponding function to write in 
* the appropriate structure
*/
  for (i=1;i<=6;i++) {
   switch (i) {
    case 1 : {
/* Hits in LGD detector*/    
     if (!LGDsave(event,nvdimlgd,nhdimlgd,nhmax)) {
      printf("ERROR-LGDsave failed\n");
      fclose(fpOut);
      return 0;
     } 
     break;          
    }
   case 2 : {
/* Hits in CPV */   
    if (!CPVsave(event,nvdimcpv,nhdimcpv,nhmax,&tindex,time_list)) {
      printf("ERROR-CPVsave failed\n");
      fclose(fpOut);
      return 0;
     } 
    break;         
   }
   case 3 : {
/* Hits in RPD */    
    if (!RPDsave(event,nvdimrpd,nhdimrpd,nhmax,&tindex,time_list)) {
      printf("ERROR-RPDsave failed\n");
      fclose(fpOut);
      return 0;
     }
    break;         
   }
   case 4 : {
/* Hits in UPV */    
    if (!UPVsave(event,nvdimupv,nhdimupv,nhmax,&tindex,time_list)) {
      printf("ERROR-UPVsave failed\n");
      fclose(fpOut);
      return 0;
     }
    break;
   }
   case 5:
/* Hits in the BSD */
     if (!BSDsave(event,nvdimbsd,nhdimbsd,nhmax,&tindex,time_list))
       {
       printf("ERROR-BSDsave failed\n");
       fclose(fpOut);
       return 0;
       }
     break;
   case 6:
/* Hits in the BGV */
     if (!BGVsave(event,nvdimbgv,nhdimbgv,nhmax,&tindex,time_list))
       {
       printf("ERROR-BGVsave failed\n");
       fclose(fpOut);
       return 0;
       }
     break;
  }              
 }
/* At this point tindex is the next available FREE entry in the time_list 
 * this gives us possibility to find the size of the time_list
*/
  size=sizeof_time_list_t(tindex);
  
/* Now add the time_list group to event */
  real_time_list=data_addGroup(event,BUFSIZE,GROUP_TIME_LIST,0,size);
  
/* Copy from time_list into real_time_list */  
  memcpy(real_time_list,time_list,size);
  
/* Do not forget to save the correct number of elements! */
  real_time_list->nelements=tindex;           

/* Now save the event */  
/* #define WRITE_ALL true */
#ifndef WRITE_ALL
  if (triggered(event))
#endif
    {
    data_clean(event);
    data_write(fileno(fpOut),event);
    }
 
/* !!!END OF THE outputEvents IS HERE!!!*/
  return 1;  
 }
