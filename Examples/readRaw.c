/*
  $Log: readRaw.c,v $
  Revision 1.4  2014/06/06 21:51:05  radphi
  * merge_clusters.c [rtj]
     - added new function to allow merging of low-mass clusters in the LGD.
  * ntuple.c [rtj]
     - implement low-mass cluster merging in the LGD if the number of photons
       is larger than 3.

  Revision 1.3  2000/01/18 03:52:40  radphi
  Changes made by jonesrt@golem
  ported Examples codes: sampleAnalysis, dumpEvent and
  paw helper routines scalers.sl, mctuple.sl, ntuple.sl to linux-Redhat (Intel)
  -rtj

 * Revision 1.2  1998/03/30  19:18:27  radphi
 * Changes made by jonesrt@
 * new implementation for March 1998 run period - R.T.Jones
 *
 * Revision 1.1  1997/02/12  15:56:18  lfcrob
 * Initial revision
 *
  */

/* Example program for how to get at raw data in an ITAPE */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <disData.h>
#include <dataIO.h>


#define BUFSIZE 100000   /* Something really big */

int main(int argc, char *argv[])
{
  FILE *fp;
  int fd;
  itape_header_t *event=NULL;
  char *inputFile=NULL;
  int *rawBuf=NULL;
  int ret;
  int nRead=0;
  int nGood=0;
  
  if(argc != 2){
    fprintf(stderr,"readRaw <inputfile>\n");    
  }
  inputFile = argv[1];
  if((fd=evt_open(inputFile,"r")) == 0){
    fprintf(stderr,"Cannot open %s\n",inputFile);
    exit(1);
  }
  
  /* Malloc the event buffer - only do this once */
  event = malloc(BUFSIZE);

  /* Note that data_read returns DATAIO_EOT at the EOD on tapes
     and at the EOF on files */

  while((ret=evt_data_read(fd,event,BUFSIZE)) != DATAIO_EOT){
    if(++nRead%1000 == 0){
      fprintf(stderr,"%d %d\r",nRead,nGood);
      fflush(stderr);
    }
    switch(ret){
    case DATAIO_OK:
      /* Now get the raw data from the ITAPE.
	 IMPORTANT: data_getGroup returns the address of the
	 group inside the ITAPE record. Thus, you NEVER EVER EVER
	 free(rawBuf) - it has not been malloc'd. The only free()
	 you would do would be on the event buffer, but their
	 is no reason to do that either, just re-use it. */	 
      
      if(((rawBuf = data_getGroup(event,GROUP_RAW,0)) == NULL) &&
         ((rawBuf = data_getGroup(event,GROUP_FROZEN,0)) == NULL)){
	/* fprintf(stderr,"No raw data!\n"); */
      }
      else{
	nGood++;
      
	/* Have fun with the raw data! The buffer rawBuf is exactly
	   the same as the buffer you would have gotten from evRead()*/
      }
      
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
  evt_close(fd);
  fprintf(stderr,"Total: %d events %d good\n",nRead,nGood);
  free(event); /* Wow, aren;t we tidy! */
  return;
}
