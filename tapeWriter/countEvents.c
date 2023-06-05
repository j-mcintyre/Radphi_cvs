/*
 * countEvents.c - count e852itape events.
 */

static const char sccsid[] = "@(#)"__FILE__"\t1.9\tCreated 6/24/97 17:27:50, \tcompiled "__DATE__;

#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <malloc.h>

#include <itypes.h>
/*#include <hdr.h>*/
#include <param.h>
#include <disData.h>
#include <tapeData.h>
#include <disIO.h>
#include <dataIO.h>
/*#include <utility.h>*/
#include <eventType.h>
#include <triggerType.h>

#define BUFSIZE 200000

int Nevents = 0;
int Tevents = 0;
float KbytesTotal = 0;
float KbytesFile = 0;

int nSync = 0;
int nPrestart = 0;
int nGo = 0;
int nEnd = 0;
int nPause = 0;
int nUnknown = 0;
int nData = 0;

int nPhysics = 0;
int nLGD = 0;
int nPed = 0;
int nBT = 0;
int nRPDMon = 0;
int nUnknownPhysics = 0;

int printAll = 0;
int keepReading = 0;
int countInc = 100;  /* Count events in multiples of 100. */
char *processName = NULL;

void ctrlCHandle(int signal);
void PrintResult(void);
int GetData(char *, FILE *finput);

PrintUsage(char *processName)
{
  fprintf(stderr,"Usage: %s [-a] file1 [file2] etc....\n\n",processName);
  fprintf(stderr,"  Options:\n");
  fprintf(stderr,"\t-a\tPrint event type sub-totals when done. (from eventType field)\n");
  fprintf(stderr,"\t-f\tContinue counting as file grows. (works like tail -f)\n");
  fprintf(stderr,"\t-h\tPrint this message.\n\n");
  exit(0);
}

main(int argc,char **argv)
{
  FILE *fp = NULL;
  int i;
  char *argptr;
  int nfile = 0,nFile,nRead;
  
  processName = argv[0];
  
  signal(SIGINT,ctrlCHandle);
  signal(SIGHUP,ctrlCHandle);
  
  fprintf(stderr,"\n");
  for (i=1; i<argc; i++) {
    argptr = argv[i];
    if (*argptr == '-') {
      argptr++;
      switch (*argptr) {
      case 'a':
	printAll = 1;
	break;
      case 'f':
	/* Note: This option may cause problems if another process is writing to
	   the file at the exact same time I am trying to read it.  When this
	   happens data_read() will usually fail because the whole event has not
	   been written out yet, which then appears to corrupt the event buffer
	   causing future data_read() calls to fail.  It would be nice if
	   data_read() could recover from these errors. */
	keepReading = 1;
	countInc = 20;  /* Decrease the event count increment. */
	break;
      case 'h':
	PrintUsage(argv[0]);
	break;
      default:
	fprintf(stderr,"Unrecognized argument: [-%s]\n\n",argptr);
	PrintUsage(argv[0]);
	break;
      }
    }
  }
  
  /* count the number of input files */
  for (nfile=0,i=1;i < argc; ++i) {
    if (*argv[i] != '-')
      nfile++;
  }
  
  if (!nfile) {
    fp = stdin;
    nfile = 1;
  }
  nFile = 1;
  nRead = 0;
  
  while (nRead < nfile) {
    if (!fp) {
      argptr = argv[nFile];
      while (*argptr == '-')
	argptr = argv[++nFile];
      nRead++;
      if (!(fp = fopen(argptr,"r")))
	fprintf(stderr,"%s: Unable to open file \'%s\': %s\n\n",argv[0],argptr,strerror(errno));
    }
    else
      argptr = "stdin";
    if (fp) {
      while (GetData(argv[0],fp)) {
	/* Do nothing - it is all done in GetData. */
      }
      fclose(fp);
      fp = NULL;
    }
    
    nFile++;
    fprintf(stderr,"\nItape records: %s has %d events (%.3f Mbytes).\n\n", argptr, Nevents, KbytesFile/1024);
    Nevents = 0;
    KbytesFile = 0;
  }
  PrintResult();
}

int GetData(char *processName, FILE *finput)
{
  static itape_header_t *buffer = NULL;
  int latch, eventType = 0;
  int ret;
  static int iloop = 3;  /* Give up after 3 EOF or errors in a row. */
  
  if (finput) {
    if (!buffer)
      buffer = (itape_header_t *)malloc(BUFSIZE);
    
    /* read from input */
    ret = data_read(fileno(finput),buffer,BUFSIZE);
    
    /* Was the read successful? */
    switch (ret) {
    default:
      fprintf(stderr,"%s: Error: data_read() failed, (ret=%d) at %d!\n",processName,ret,Nevents);
      iloop--;
      if (keepReading)
	sleep(2);
      break;
    case DATAIO_EOT:
      if (keepReading) {
	sleep(2);  /* So I don't waste CPU time. */
      } else {
	fprintf(stderr,"%s: EOF at %d (%.3f Mbytes)...\n", processName, Nevents, KbytesTotal/1024);
	return(0);
      }
      break;
    case DATAIO_EOF:
      if (keepReading) {
	sleep(2);  /* So I don't waste CPU time. */
      } else {
	fprintf(stderr,"%s: EOT at %d (%.3f Mbytes)...\n", processName, Nevents, KbytesTotal/1024);
	iloop--;
      }
      break;
    case DATAIO_OK:
      iloop = 3;
      
      /* Is this real data? */
      switch (buffer->type) {
      default:
	fprintf(stderr,"%s: Unknown data type: %d.\n",processName,buffer->type);
	break;
      case TYPE_TAPEHEADER:
	{
	  tapeHeader_t *hdr = (void*)buffer;
	  fprintf(stderr,"\nReading tape number %04d, label written on %s\n",hdr->tapeSeqNumber,ctime((void*)&hdr->startTime));
	}
      break;
      case TYPE_ITAPE:
	Nevents++;  /* Only want to count real data */
	Tevents++;
	KbytesFile += ((float)buffer->length + 4) / 1024;  /* extra 4 bytes of CRC per event */
	KbytesTotal += ((float)buffer->length + 4) / 1024;  /* extra 4 bytes of CRC per event */
	if (!(Tevents % countInc)) {
	  fprintf(stderr,"%d \r",Tevents);
	  fflush(stderr);
	}
	if (printAll) {
	  eventType = buffer->eventType;
	  switch(eventType){
	  case EV_DATA:
	    nData++;
	    switch(buffer->trigger){
	    case TRIG_DATA:
	      nPhysics++;
	      break;	      
	    case TRIG_LGDMON:
	      nLGD++;
	      break;
	    case TRIG_PED:
	      nPed++;
	      break;
	    case TRIG_BT:
	      nBT++;
	      break;
	    case TRIG_RPDMON:
	      nRPDMon++;
	      break;
	    }
	    break;
	  case EV_SYNC:
	    nSync++;
	    break;
	  case EV_PRESTART:
	    nPrestart++;
	    break;
	  case EV_GO:
	    nGo++;
	    break;
	  case EV_PAUSE:
	    nPause++;
	    break;
	  case EV_END:
	    nEnd++;
	    break;
	  case EV_UNKNOWN:
	    nUnknown++;
	    break;
	  }
	}
	break;
      }
      break;
    }
  }
  return(iloop);
}

void ctrlCHandle(int signal)
{
  /* signal(SIGINT,ctrlCHandle);
  signal(SIGHUP,ctrlCHandle); */
  fprintf(stderr,"\n%s: ",processName);
  PrintResult();
  exit(1);
}

void PrintResult(void)
{
  if(printAll) {
    fprintf(stderr,"Sub-Totals: %d Control, %d Data Triggers.\n\n",nSync+nGo+nEnd+nPause+nPrestart,nData);
    fprintf(stderr,"\tControl Triggers\tData Triggers\n");
    fprintf(stderr,"\t-------------------------------------------------\n");
    fprintf(stderr,"\tPrestart\t %d\tPhysics\t\t%d\n",nPrestart,nPhysics);
    fprintf(stderr,"\tGo\t\t %d\tLGD\t\t%d\n",nGo,nLGD);
    fprintf(stderr,"\tPause\t\t %d\tPedestal\t%d\n",nPause,nPed);
    fprintf(stderr,"\tSync\t\t %d\tBase Test\t%d\n",nSync,nBT);
    fprintf(stderr,"\tEnd\t\t %d\tRPD Mon\t\t%d\n",nEnd,nRPDMon);
  }
}

/* end file */
