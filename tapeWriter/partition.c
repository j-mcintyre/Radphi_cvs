/*
 * partition
 *
 * A standalone program to input a stream of data on standard input and
 * split the content output into multiple output streams.  The input and
 * output streams must be in itape format.  The data from the input stream
 * is divided into TRIG_DATA records and the rest, called special records.
 * Special records are copied verbatim into partition 0.  TRIG_DATA records
 * are copied into a series of partitions starting at partition 1.  Each
 * partition can hold up to partitionSize records, where partitionSize can
 * be set on the command line or left at the default value set below.  The
 * output files are named <basename>p0.itape ... <basename>pN.itape where
 * N is determined by the partitionSize and the length of the input stream.
 *
 * Richard Jones, July 10 2000
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <eventType.h>
#include <triggerType.h>
#include <disData.h>
#include <dataIO.h>
#include <itypes.h>
#include <tapeData.h>
#include <math.h>
#include <unpackData.h>

#define BUFSIZE 100000

itape_header_t *event;
int partitionSize = 100000;

void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  partition <options> basename\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-h\t\tPrint this message\n");
  fprintf(stderr,"\t-p#\t\tStore a most # records per partition.\n");
}

int main(int argc, char *argv[])

{
  int iarg;
  char *argptr;
  int fdOut[2];
  int nEvents;
  int partition;
  char basename[80];
  char fullname[100];

  event = malloc(BUFSIZE);
  if (event == 0) {
    fprintf(stderr,"Failed to allocate event buffer of %d bytes\n",BUFSIZE);
    exit(1);
  }

  tape_setStream(stdin);

  for (iarg=1;iarg<argc;iarg++) {
    if (*argv[iarg] == '-') {
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'p':
	partitionSize = atoi(++argptr);
	break;
      case 'h':
	Usage();
	exit(1);
      default:
	fprintf(stderr,"Unknown argument: %s\n",argv[iarg]);
	Usage();
	exit(1);
      }
    }
    else {
      strcpy(basename,argv[iarg]);
    }
  }
  
  nEvents = partition = 0;

  sprintf(fullname,"%sp%d.itape",basename,partition++);
  if ((fdOut[0] = evt_open(fullname,"w")) == 0) {
    fprintf(stderr,"Failed to open %s for writing\n",fullname);
    exit(1);
  }

  while (data_read(fileno(stdin),event,BUFSIZE) == DATAIO_OK) {
    if (event->eventType == EV_DATA && event->trigger == TRIG_DATA ||
        event->eventType == EV_MCGEN) {
      if (partitionSize > 0) {
        if ((++nEvents % partitionSize) == 1) {
          if (nEvents > 1) {
            evt_close(fdOut[1]);
          }
          sprintf(fullname,"%sp%d.itape",basename,partition++);
          if ((fdOut[1] = evt_open(fullname,"w")) == 0) {
            fprintf(stderr,"Failed to open %s for writing\n",fullname);
            exit(1);
          }
          nEvents = 1;
        }
        data_write(fdOut[1],event);
      }
    }
    else {
      data_write(fdOut[0],event);
    }
  }
  fflush(NULL);
  return 0;
}
