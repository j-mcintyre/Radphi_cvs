#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <dataIO.h>
#include <disData.h>

#define BUFSIZE 100000

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Dummy call to load the module
*/
void fileIO(void) {}    /* c entry point */
void fileio_(void) {}   /* fortran entry */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*  fileOpen returns file descriptor on success, otherwise 0.
*/
int fileOpen(char *fname, char *pmode){
  int fd=0;
  char *fin;
  if ((fin=(char *)index(fname,' ')) != NULL)
    fin[0]=0;
  if ((fd=evt_open(fname,pmode)) == 0) {
    fprintf(stderr,"ERROR - failed to open input event file \"%s\"\n",fname);
    return -1;
  }
  return fd;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*  fileGetEvent returns eventNo if DATAIO_OK, -1 if error,
*  and 0 for endOfTape.
*/
int fileGetEvent(int fd, itape_header_t *event) {
   int res = evt_data_read(fd,event,BUFSIZE);
   switch(res) {
   case DATAIO_ERROR:
     fprintf(stderr,"fileGetEvent: data_read() returned DATAIO_ERROR!\n");
     return -2;
   case DATAIO_BADCRC:
     fprintf(stderr,"fileGetEvent: data_read() returned DATAIO_BADCRC!\n");
     return -3; 
   case DATAIO_EOT:
     return 0;
   case DATAIO_EOF:
     return fileGetEvent(fd,event);
   case DATAIO_OK:
     return event->eventNo;
   default:
     fprintf(stderr,"fileGetEvent: unknown code returned by data_read()!\n");
     return -9;
   } 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* fileSkipEvent skip desired number of events (records)
*/
int fileSkipEvents(int fd, int count) {
  
  itape_header_t *event=NULL; 
  int rec=0;
  int ret;

  event = malloc(BUFSIZE);
  if (event == NULL) {
    fprintf(stderr,"ERROR - failed to allocate aux event buffer!\n");
    return -9;
  }
  
  fprintf(stderr," Skipping %d events ... \n",count);
  while ( (rec < count) && (ret=fileGetEvent(fd,event) > 0) ) {rec++;}
  if ( rec < count ) {
    fprintf(stderr,"Warning - skipped %d instead of %d events\n",rec,count);
  }

  free(event);
  return rec;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*  Close file for a given file descriptor.
*/
void fileClose(int fd){
  if ( fd <= 0) {
    fprintf(stderr,"WARNINIG - Invalid file descriptor! \n");
    return;
  }
  evt_close(fd);
  return;
}

int filetest_(char *infile) {

  int fd=0;
  int ret=0;
  int record=0;                 /* current record in the event bufer */
  int records=0;                /* should be = record for the same file */
  char *pmode = "r";
  int skip = 15;
  int start = 5;
  
  itape_header_t *event=NULL;

/* allocate event */
  if (event == NULL) {
    event = malloc(BUFSIZE);
    if (event == NULL) {
      fprintf(stderr,"ERROR - failed to allocate event buffer!\n");
      return -9;
    }
  }

/* * * * Open file */
  if ( (fd=fileOpen(infile,pmode)) < 0 ) {
    fprintf(stderr,"Exit after openItape failed on %s \n",infile);
    return -1;
  }

/* * *  skip events  */
  if ( (records=fileSkipEvents(fd,skip)) != skip  ) {
    fprintf(stderr,"WARNING - skipped less than %d,",skip);
    fprintf(stderr,"probably passed EOF\n");
  }
  fprintf(stderr,"* * * * * after skippping %d events ",records);
  fprintf(stderr,"ready to read event %d\n",(records+record+1));

  if ( (record=fileGetEvent(fd,event)) < 0)  {
    fprintf(stderr,"ERROR - exit after reading event \n");
    return record;
  }
  fprintf(stderr,"* * * * *  read event %d  \n",record);
  fprintf(stderr,"* * * * *  run number was %d \n",event->runNo);

  fileClose(fd);
  return records;
}
