/* fetches MonteCarlo event stubs from an event generator file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <disData.h>
#include <dataIO.h>
}

#include <rootclasses.h>

#define BUFSIZE 100000   /* Something really big */

/* process the next event from the file */
int processEvent(char *filename)
{
  static FILE *fp=NULL;
  static itape_header_t *event=NULL;
  int ret;
  
  if ((fp == NULL) && ((fp=fopen(filename,"r")) == NULL)) {
    fprintf(stderr,"WARNING - failed to open input events file %s\n",filename);
    return 1;
  }
  
  /* Malloc the event buffer - only do this once */
  if ((event == NULL) && ((event=(itape_header_t *)malloc(BUFSIZE)) == NULL)) {
    fprintf(stderr,"ERROR - failed to allocating input event buffer!\n");
    fclose(fp);
    return 1;
  }

  TMCevent mc_event(event,BUFSIZE);
  TLGDhits lgd_hits(event,BUFSIZE);
  TCPVhits cpv_hits(event,BUFSIZE);
  TUPVhits upv_hits(event,BUFSIZE);
  TBSDhits bsd_hits(event,BUFSIZE);
  TBGVhits bgv_hits(event,BUFSIZE);
  TTimeList time_list(event,BUFSIZE);

  /* Note that data_read returns DATAIO_EOT at the EOD on tapes
     and at the EOF on files */

  if ((ret=data_read(fileno(fp),event,BUFSIZE)) == DATAIO_OK) {
    printf("===>New event, run %d, event %d\n",event->runNo,event->eventNo);
    
    if (mc_event.GetGroup() == NULL) {
      fprintf(stderr,"Warning: mc_event group missing event %d\n",event->eventNo);
    }
    else {
      mc_event.Print();
    }
    
    if (lgd_hits.GetGroup() == NULL) {
      fprintf(stderr,"ERROR - lgd_hits group missing from input data\n");
    }  
    else {
      lgd_hits.Print();
    }
 
    if (cpv_hits.GetGroup() == NULL) {
      fprintf(stderr,"ERROR - cpv_hits group missing from input data\n");
    } 
    else {
      cpv_hits.Print();
    }

    if (upv_hits.GetGroup() == NULL) {
      fprintf(stderr,"ERROR - upv_hits group missing from input data\n");
    }
    else {
      upv_hits.Print();
    } 
  
    if (bsd_hits.GetGroup() == NULL) {
      fprintf(stderr,"ERROR - bsd_hits group missing from input data\n");
    }
    else {
      bsd_hits.Print();
    }

    if (bgv_hits.GetGroup() == NULL) {
      fprintf(stderr,"ERROR - bgv_group missing from input data\n");
    }
    else {
      bgv_hits.Print();
    }

    if (time_list.GetGroup() == NULL) {
      fprintf(stderr,"ERROR - time list group missing from input data\n");
    } 
    else {
      time_list.Print();
    }
  }
  else if (ret!=DATAIO_EOT) {
    fprintf(stderr,"ERROR - Couldn't read event group from file %s\n",filename);
    fclose(fp);
    return 1;
  }
  else {
   fclose(fp);
   return 2;
  }

  return 0;
}

/* main program */
int main(int argc, char *argv[])
{
  int i=1;
  char filename[120];
  if ((argc>2)||((argc==2)&&(argv[1][0]=='-'))) {
    printf("usage: %s [filename]\n",argv[0]);
    exit(1);
  }
  if (argc==1)
    strcpy(filename,"simData.itape");
  else
    strcpy(filename,argv[1]);
  while (!processEvent(filename)) {
    printf("==============================================================\n");
    printf("This was event number %i\n",i++);
    printf("==============================================================\n");
  }
  return 0;
} 
