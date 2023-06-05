#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include "runs.h"
#include "errors.h"
#include "global.h"

#include <cernlib.h>

#include <disData.h>
#include <dataIO.h>
#include <ntypes.h>
#include <itypes.h>
#include <tapeData.h>
#include <umap.h>
#include <lgdGeom.h>
#include <pedestals.h>
#include <eventType.h>
#include <triggerType.h>
#include <disIO.h>
#include <unpackData.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <clib.h>
#include <lgdCal.h>
#include <calibration.h>
#include <makeHits.h>
#include <libdecode.h>

#define PREFILL 10

static int getAnEvent(itape_header_t *event);
static int dispatcherReconnect(const char*host,int pipelinePrefill);

void setupHbook(void);

static int in_tape_fd, out_tape_fd;      /* file descriptors */
static FILE *in_tape, *out_tape;         /* file handles for I/O */
/* char *histo_file_name = "radphi.hbook"; */

static float *ped;
static int *thresh;
static float *cc;
static int nChannels,nRows,nCols;
static float baseCC=0.00175;

itape_header_t *tape = NULL;     /* tape record we will write to  */

void analysis() {
  static int first=1;
  int size;                        /* used to calculate size of groups before adding to tape */

  itape_header_t *event = NULL;    /* tape record coming from input */
  int *rawbuf=NULL;                /* pointer into raw data group (CEF) */
  int *group;                      /* pointer returned from data_addGroup */

  int length;                 
  int32 run;                       /* run number from list of runs to process (runs.h) */
  int ROCnum, slot;                /* current ROC number and current geographical address */
  int32 upper, lower;              /* upper and lower bounds of event numbers in run (runs.h) */
  int32 db_run = -1;               /* run number of currently loaded hardware database */

  int ERR = 0;                     /* flag that gets set if we have errors in reading raw data */

  time_t begin = time(NULL);       /* mark entry to processing */
  int EV_DATA_count = 0;           /* used to count number of successful EV_DATA events */
  int TRIG_count[6];               /* counter for physics trigger types  */
  int i;


  lgd_hits_t *lgdHits=NULL,*clusterHits=NULL;  
  lgd_clusters_t *clusters=NULL, *tmp=NULL;


  for (i=0; i<6; i++) TRIG_count[i]=0;

  if(first){
    unpackSetLogFile(logfi);
    setupHbook();
    first=0;
  }

  /* get all the data tapes ready */

  if(dispatcher){
    if(dispatcherReconnect(dispatcher,PREFILL)){
      fprintf(stderr,"Failed to connect to dispatcher!!\n");
      fprintf(logfi,"Failed to connect to dispatcher!!\n");
      exit(1);
    }      
  }
  else{
    if (input_tape == NULL) {
      in_tape_fd = 0;
      fprintf(logfi, "Input coming from stdin\n");
    }
    else {
      if ((in_tape_fd = evt_open(input_tape, "r")) == 0)
	file_error(input_tape);
      fprintf(logfi, "Input file: %s\n", input_tape);
    }
    
    if (output_tape == NULL) {
      out_tape_fd = 1;
      fprintf(logfi, "Output going to stdout\n"); 
    }  
    else if (strcmp(output_tape,"/dev/null") == 0) {
      out_tape_fd = 0;
      fprintf(logfi, "Output turned off\n"); 
    }
    else {
      if ((out_tape_fd = evt_open(output_tape, "w")) == 0)
	file_error(output_tape);
      fprintf(logfi, "Output file: %s\n", output_tape);
    }
  }
    
  /* allocate some ITAPE records */

  if ((event = (itape_header_t *)malloc(100000)) == NULL )
    malloc_error(__FILE__, __LINE__);
  if ((tape = (itape_header_t *)malloc(100000)) == NULL )
    malloc_error(__FILE__, __LINE__);
 


  while(get_next_node(&run, &lower, &upper)) {


    while(getAnEvent(event) == 0){

      /* 
       *  we want to know what kind of event we have 
       *  for now, all of the events are CODA events
       *  just mention that we found them
       */

      switch(event->eventType) {
      case EV_SYNC:
	fprintf(logfi, "  Found CODA Sync event\n");
	break;
      case EV_PRESTART:
	fprintf(logfi, "  Found CODA Prestart event\n");
	break;
      case EV_GO:
	fprintf(logfi, "  Found CODA Go event\n");
	break;
      case EV_PAUSE:
	fprintf(logfi, "  Found CODA Pause event\n");
	break;
      case EV_END:
	fprintf(logfi, "  Found CODA End event\n");
	break;
      case EV_UNKNOWN:
	fprintf(logfi, "  Found an unknown CODA event type\n");
	break;
      case EV_DATA: break;   /* go ahead and process it */
      default: break;        /* this event is NOT raw data */
      }

      /* check to see if we are in the bounds of what we want to process */

      lower=0;
      if (event->eventNo < lower) continue;
      else if (event->eventNo > upper) {
	if (!get_next_node(&run, &lower, &upper)) break;
	else if (event->eventNo < lower) continue;
	else if (event->eventNo > upper) break;
      }

      /* lets process this event shall we? */
      /* new record on tape */
      if (data_newItape(tape)) {
	fprintf(logfi, "ERROR: failed to allocate new Itape record\n");
	fprintf(stderr, "ERROR: failed to allocate new Itape record\n");
	fprintf(stderr, ">>>>>> Current event: %d, run: %d\n", event->eventNo, event->runNo);
	clear_runs();
	break;
      }

      if (out_tape_fd) {
        tape->type = TYPE_ITAPE;
        tape->eventNo = event->eventNo;
        tape->time = time(NULL);    /* mark this event with the current time */
        tape->runNo = event->runNo;
        tape->eventType = event->eventType;   /* for now, just store EV_DATA */
        tape->trigger = event->trigger;
      
        /* Here we will output the GROUP_RAW, if it exists  */

        if((rawbuf=data_getGroup(event, GROUP_FROZEN, 0)) != NULL) {
	  length = *rawbuf;
	  size = (length + 1)*sizeof(uint32);
	  group = data_addGroup(tape, 100000, GROUP_FROZEN, 0, size);
	  memcpy(group, rawbuf, size);
        }
        if((rawbuf= data_getGroup(event, GROUP_RAW, 0)) != NULL) {
	  length = *rawbuf;
	  size = (length + 1)*sizeof(uint32);
	  group = data_addGroup(tape, 100000, GROUP_RAW, 0, size);
	  memcpy(group, rawbuf, size);
        }

        /*  Here we will output Monte Carlo data, if it exists  */

        if (data_getGroup(event, GROUP_MC_EVENT, 0) != NULL) {
          int vertex;
          for (vertex=0;
               rawbuf=data_getGroup(event,GROUP_MC_EVENT,vertex);
               vertex++) {
	    length = *rawbuf;
	    size = (length + 1)*sizeof(uint32);
	    group = data_addGroup(tape, 100000, GROUP_MC_EVENT, vertex, size);
	    memcpy(group, rawbuf, size);
          }
	  monte_carlo = 1;
	  event->eventType = EV_DATA; /* force all MC events as real data */
        }
      }

      /* Actually unpack the event */

      if (event->eventType == EV_DATA) {
	if (event->runNo != db_run) {
	  db_run = event->runNo;
	  if(setup_makeHits(event)){
	    fprintf(stderr,"Problem with setup_makeHits; exiting\n");
	    fprintf(logfi,"Problem with setup_makeHits; exiting\n");
	    return;
	  }      
	}
	if(monte_carlo || !unpackEvent(event,BUFSIZE)){
          tdc_values_t *bsd_tdcs = data_getGroup(event,GROUP_BSD_TDCS,0);
	  tdc_values_t *cpv_tdcs = data_getGroup(event,GROUP_CPV_TDCS,0);
	  tdc_values_t *upv_tdcs = data_getGroup(event,GROUP_UPV_TDCS,0);
	  tdc_values_t *bgv_tdcs = data_getGroup(event,GROUP_BGV_TDCS,0);
	  tdc_values_t *tagger_tdcs = data_getGroup(event,GROUP_TAGGER_TDCS,0);
          if (!bsd_tdcs && !tagger_tdcs && !bgv_tdcs
                        && !cpv_tdcs && !upv_tdcs) {
            /* ignore events without tdc info (1 every 10000) */
            continue;
          }
         /* don't process unless physics trigger or Monte Carlo */
	  if(event->trigger ==2 || monte_carlo){
	    makeHits(event);
            if (process(event)) {
	      if(out_tape_fd){
	        data_write(out_tape_fd, tape);
             /*	data_flush(out_tape_fd);  not on every event, please! -rtj- */
	      }
	      else if(dispatcher){
	        data_addCRC(event);
	        disIO_writeCOOKED(event,event->length+4);
	      }
            }
          }

	  EV_DATA_count++;
          {
             double count;
             for (count = EV_DATA_count; count > 9; count /= 10) {}
             if (count == (int)count)
	        fprintf(stderr," processing event  %d\r",EV_DATA_count);
	  }
          
	  {
	    int i=0;
	    int tmp=event->trigger;
	    while(tmp >> i != 1)
	      i++;
	    
	    TRIG_count[i]++;
	  }
	}
      }    /* end EV_DATA */
	
      else if (event->eventType != EV_DATA &&
               event->eventType != EV_PRESTART &&
               event->eventType != EV_GO &&
               event->eventType != EV_END) {
	  /*
	   *  In CODA, it is possible to define 16 different event types
	   *  These types are encoded in the event headers
	   *  These different event types can then be treated differently in 
	   *  offline data analysis.
	   *  Currently, there is only one event type, anything else is
	   *  treated as an error
	   */
	  if(dispatcher){
	      data_addCRC(event);
	      disIO_writeCOOKED(event,event->length+4);
	  }
	  fprintf(stderr, "Unknown raw event type (0x%x)!\n",event->eventType);
	  fprintf(logfi, "ERROR: Encountered a raw event type not yet defined\n");
	  fprintf(logfi, "       Run: %d, Event:%d\n", event->runNo, event->eventNo);
      }
    }     /* end while data read */
  } /* end while get next node */

  /* we're done, cleanup and close down */

  free(event);
  free(tape);
  if(dispatcher){
    disIO_disconnect();
  }
  else{
    if (in_tape_fd) evt_close(in_tape_fd);
    if (out_tape_fd) evt_close(out_tape_fd);
  }
  fprintf(logfi, "****** Ending data analysis. ********\n\nTime taken for analysis: %10.2f seconds\n", difftime(time(NULL), begin));
  fprintf(logfi, "  Number of correctly parsed events of EV_DATA: %d from run %d\n", EV_DATA_count,event->runNo);
  fprintf(stderr, "  Number of correctly parsed events of EV_DATA: %d from run %d\n", EV_DATA_count,event->runNo);

  /* modify the below for data starting in summer 1999 when we went to 5 trigger types */
  /*   for (i=1; i<6; i++){
      fprintf(logfi, "  Number of correctly parsed events of trigger %d : %d from run %d\n", i,TRIG_count[i],event->runNo);
      fprintf(stderr, "  Number of correctly parsed events of trigger %d : %d from run %d\n", i,TRIG_count[i],event->runNo);

   }
   */

  if (event->runNo > 6000){
      fprintf(logfi, "  Number of correctly parsed triggers from run %d:\n", event->runNo);
      fprintf(stderr, "  Number of correctly parsed triggers from run %d:\n", event->runNo);

      fprintf(stderr, "      Physics triggers   : %d\n      LGD laser Monitor  : %d\n      Barrel Laser       : %d\n      pedestal           : %d\n      basetest & scalers : %d\n ",  TRIG_count[1],TRIG_count[2],TRIG_count[3],TRIG_count[4],TRIG_count[5]);

      fprintf(logfi, "      Physics triggers   : %d\n      LGD laser Monitor  : %d\n      Barrel Laser       : %d\n      pedestal           : %d\n      basetest & scalers : %d\n ",  TRIG_count[1],TRIG_count[2],TRIG_count[3],TRIG_count[4],TRIG_count[5]);
  }
  else{
      fprintf(logfi, "  Number of correctly parsed triggers from run %d:\n", event->runNo);
      fprintf(stderr, "  Number of correctly parsed triggers from run %d:\n", event->runNo);

      fprintf(stderr, "      Physics triggers   : %d\n      LGD laser Monitor  : %d\n      Pedestal           : %d\n      Basetest & scalers : %d\n  ",  TRIG_count[1],TRIG_count[2],TRIG_count[3],TRIG_count[4]);
      fprintf(logfi, "      Physics triggers   : %d\n      LGD laser Monitor  : %d\n      Pedestal           : %d\n      Basetest & scalers : %d\n  ",  TRIG_count[1],TRIG_count[2],TRIG_count[3],TRIG_count[4]);
  }
  return;
}

static int getAnEvent(itape_header_t *event)
{
  int index;
  int ret=DISIO_EOF;
  int raw_size;
  char *raw_data=NULL;
  tapeHeader_t *hdr=NULL;
  if(dispatcher){
    disIO_command("REQUEST_DATA");
    ret =  disIO_readRAW(event,100000,&raw_size,1);    
    if(ret<0){ 
      fprintf(stderr,"from_dispatcher: Error: disIO_readRAW() returned %d\n",ret);
      exit(1);
    }
    raw_data = (char *)event;
    switch(ret){
    default:
    case DISIO_COMMAND:      
      fprintf(stderr,"Command from Dispatcher: %s (%d)\n",raw_data,raw_size);
      {
	const char *ptr = strtok(raw_data,":");
	if(!strcmp(ptr,"QUIT")){
	  return(0);
	}
	if(!strcmp(ptr,"BEGINTAPE")){
	  for(index=0;index<PREFILL;index++){
	    disIO_command("REQUEST_DATA");
	  }
	}
	if(!strcmp(ptr,"ENDTAPE")){
	  disIO_command("FINISHED");
	  return(getAnEvent(event));
	}
	return(getAnEvent(event));
      }
    case DISIO_OK:
      break;
    }
  }
  else{
    ret = evt_data_read(in_tape_fd, event, 100000);
    if(ret != DATAIO_OK)
      return(1);
  }
  /*  here deal with itape headers... */
    /* Got an event */
    switch(event->type){
    case TYPE_TAPEHEADER:
      hdr = (void *) event;
      return(getAnEvent(event));
    default:
      return(0);
    }
}

static int dispatcherReconnect(const char*host,int pipelinePrefill)
{
  int i;
  int retry = 2;
  int sleepTime = 30; /* seconds between retries */

  if (disIO_socket >= 0){
    disIO_command("QUIT");
    disIO_disconnect();
  }
	  
  while(retry--){
    int ret = disIO_connect(host,0);
    if (ret==0){
      fprintf(stderr,"Connected to the Dispatcher at %s\n",host);	
      for (i=0; i<pipelinePrefill; i++){
	disIO_command("REQUEST_DATA");
      }	
      return 0;
    }
    fprintf(stderr,"Cannot connect to the Dispatcher at %s, waiting %d sec before retry...\n",host,sleepTime);
    sleep(sleepTime);
  }
  fprintf(stderr,"Cannot connect to the Dispatcher at %s, giving up.\n",host);
  return -1;
}
