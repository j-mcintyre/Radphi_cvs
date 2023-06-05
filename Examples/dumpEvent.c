/* dumps the contents of an ITAPE file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <disData.h>
#include <dataIO.h>


#define BUFSIZE 100000   /* Something really big */

int main(int argc, char *argv[])
{
  int fp;
  itape_header_t *event=NULL;
  char *inputFile=NULL;
  int ngroups=0;
  uint32 *types;
  uint32 *subtypes;
  size_t *sizes;
  int ret;
  int i;
  
  if(argc != 2){
    fprintf(stderr,"Usage: dumpEvent <inputfile>\n");    
    exit(1);
  }
  inputFile = argv[1];
  if((fp=evt_open(inputFile,"r")) == 0){
    fprintf(stderr,"Cannot open %s\n",inputFile);
    exit(1);
  }
  
  /* Malloc the event buffer - only do this once */
  event = malloc(BUFSIZE);

  /* Note that data_read returns DATAIO_EOT at the EOD on tapes
     and at the EOF on files */

  while((ret=evt_data_read(fp,event,BUFSIZE)) != DATAIO_EOT){
    switch(ret){
    case DATAIO_OK:
      printf("\n>>--->New ITAPE event\n");
      printf("length %10d type %10d transComputerCode %10d ",
              event->length,event->type,event->transComputerCode);
      printf("run %5d spill %10d event %10d\n",
              event->runNo,event->spillNo,event->eventNo);
      printf("eventType %10d trigger %10d time %10d latch %10d\n",
              event->eventType,event->trigger,event->time,event->latch);
      if (data_listGroups(event,&ngroups,&types,&subtypes,&sizes)){
        fprintf(stderr,"dumpEvent - data_listGroups error, aborting!\n");
        exit(1);
      }
      for (i=0;i<ngroups;i++)
        printf("group %5d type %10d subtype %10d length %10d\n",
                i+1,*(types++),*(subtypes++),*(sizes++));
      break;
    case DATAIO_EOF:
      /* This is just a file marker on a tape - read past it */
      fprintf(stderr,"EOF marker!\n");
      break;
    case DATAIO_ERROR:
      /* Bad, try to continue */
      fprintf(stderr,"Error reading data!\n");
      break;
    case DATAIO_BADCRC:
      /* Bad CRC, event is hosed */
      fprintf(stderr,"Error reading data (bad CRC)!\n");
      break;
    default:
      /* Something REALLY bad - data_read completely failed */
      fprintf(stderr,"Unknown return code form data_read : %d\n",ret);
      break;
    }
  }
  free(event); /* Wow, aren;t we tidy! */
  return;
}
