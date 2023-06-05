/*
 * tapeIO.h
 *
*/

#ifndef TAPEIOH_INCLUDED
#define TAPEIOH_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <ntypes.h>

#define TAPEIO_OK       0
#define TAPEIO_ERROR  (-1)
#define TAPEIO_EOT    (-2)   /* end of data (End Of Tape) */

#define TAPEIO_STREAM   1    /* file descriptor is not a tape device      */
#define TAPEIO_FIXED    2    /* it's a tape device in fixed block mode    */
#define TAPEIO_VARIABLE 3    /* it's a tape device in variable block mode */

int tape_devType(int fd); /* returns one of the TAPEIO_xxxx or (-1) if error */

FILE *tape_getStream(int fd);
int tape_setStream(FILE* fp);
int tape_vwrite(int p);
int tape_vread(int p);
int tape_flush(int fd);
int tape_ff(int fd,off_t skipbytes);
int tape_fwrite_fillBuf(int p,const void *ptr,int len);
int tape_write(int fd,const void*buffer,int len);

int tape_getBlockSize(int fd,uint32*minSize,uint32*maxSize,uint32*recSize);
int tape_getState(int fd);              /* returns bits defined in mtio.h: MT_BOT, etc... */
uint32 tape_getPosition(int fd); /* returns current tape block number (does not work for jag tapes) */

int tape_setBlockSize(int fd,uint32 newBlockSize);

int tape_rewind(int fd);
int tape_unload(int fd);
int tape_findFM(int fd,int count);  /* count > 0 seeks forward, count < 0 seeks backwards */
int tape_writeFM(int fd);

#endif
/* end file */
