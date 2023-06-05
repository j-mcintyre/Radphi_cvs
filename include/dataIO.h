/*
 * dataIO.h
 *
*/

#include <stdio.h>
#include <stdlib.h>

#ifndef _DATAIO_INCLUDED
#define _DATAIO_INCLUDED

#define DATAIO_EOF      0               /* file mark or end of tape              */
#define DATAIO_OK       1               /* no errors                             */
#define DATAIO_EOT      2               /* end of tape                           */
#define DATAIO_ERROR  (-1)              /* read/write error (message->stderr)    */
#define DATAIO_BADCRC (-2)              /* CRC mismatch error: data is damaged   */

/* 
   The evt_open and evt_data_read commands are provide a uniform i/o
   interface to raw CODA and itape event files.  Using them instead of
   the normal fopen and data_read sequence avoids the need to run the
   wrapFile command on every file before looking at its contents.  A
   file opened with the evt_fopen command should be closed with evt_close.
*/
int evt_open(char *filename,char *mode);

int evt_data_read(int fd,void *buffer,int bufsize);

int evt_close(int fd);

int data_write(int fd,const void* event);         /* automatically regenerates the CRC word        */

int data_read(int fd,void* buffer,int bufsize);
int data_read_alloc(int fd,void **buffer);         /* memory pointed by *buffer should be free()ed  */

int data_flush(int fd);   /* Flush internal buffers: should always be called before a close()       */

int data_writeFM(int fd);                          /* Write a file mark on tape                     */

int data_findFM(int fd,int count); /* find a File Mark. count>0 seeks forward, count<0 seeks backwards */

/*
 * data_SeekTape: Position the data tape at the requested run number.
 *
 * Input arguments:     fd - file descriptor of the tape device (as returned by 'open' or 'fileno')
 *                   runNo - requested run number
 *
 * Return value:    0  if run was found.
 *                (-1) if run was not found
 *
 * Output arguments:  tapeNo - if a tape label is encountered while
 *                             seeking the requested run,
 *                             'tapeNo' will be set to the run number from it.
 *                             Otherwise the value is left unchanged.
*/

int data_SeekTape(int fd,int runNo,int *tapeNo);

/*
 * data_setSieve: declare an output file as a sieve record of an
 *                existing input file
 *
 * input arguments:
 * 	fp 	 - stream pointer to open sieve file
 * 	fdin	 - descriptor of open input file (see filename below)
 * 		   if writing to fp, or 0 if reading from fp
 *	filename - pathname of associated input file (see fdin above)
 *	           if writing to fp, or 0 if reading from fp
 *
 * output arguments: none
 *
 * return value:  0  if success
 *               -1  if failure
 */

int data_setSieve(FILE *fp, int  fdin, char* filename);

/*
 * data_ff: fast-forward over a fixed number of bytes on input stream;
 *          duplicates functionality of data_findFM for streams.
 *
 * input arguments:
 * 	fd	  - file descriptor (must be open to a stream)
 * 	skipbytes - byte count to be skipped
 *
 * output arguments: none
 *
 * return value:  0  if success
 *               -1  if failure
 */

int data_ff(int fd, off_t skipbytes);

#endif
/* endfile */
