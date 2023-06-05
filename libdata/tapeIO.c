/*
 * tapeIO.c
 *
 * Modification history
 *
 * July 30, 2004 - Richard Jones
 * - fixed a bug in the stream_reader where it overcounted the nbytes in
 *   dummy buffers while doing a fast-forward operation.  This happened
 *   because the bufdes->filemark could be increased by tape_ff in the
 *   main thread while the stream_reader was doing a fseek to the original
 *   offset. One should always make a local copy of mutable globals before
 *   doing non-atomic actions based on their values!  While hunting that
 *   one, added a journaling facility to help debugging of multi-threaded
 *   i/o, improved semantics of the buffer interlocks between the threads,
 *   added a field in each buffer showing the stream offset to its data. 
 *
 * February 17, 2003 - Richard Jones
 * - added tape_ff to provide fast-forward capability on input streams.
 *   This functionality is needed to implement data sieves in the dataIO
 *   package.  To implement tape_ff with multibuffer readahead I had to
 *   introduce new members bytecount and filemark to the tape_info struct.
 *   The bytecount is a 64-bit running count of all bytes transfered on
 *   a descriptor.  It may be of more general use at some point, so I
 *   implemented the counter on all units, both input and output.  
 *
 * October 19, 2001 - Richard Jones
 * - fixed a bug in the multi-buffer io flow diagram using the condition
 *   signal facility of the pthreads library.  In the old code the io
 *   thread was blocked from overtaking the user thread by a mutex held
 *   by the user thread on the buffer being processed.  But if the io
 *   thread got blocked after it released a mutex and before it acquired
 *   the next one then the user thread could buzz right past the io thread
 *   and process buffers containing no (or stale) data.  This was not
 *   infrequent.  Now the user thread explicitly checks each buffer status
 *   to see that it is filled before using it.  It uses pthread_cond_wait
 *   to block if the buffer it acquired is not yet filled.
 * 
 * May 10, 2001 - Richard Jones
 * - added overlapping of reads from a stream for high-efficiency i/o.
 *   This requires the pthreads library.  It has been tested under
 *   gcc version egcs-2.91.66 19990314/Linux (egcs-1.1.2 release)
 *   and is enabled by defining the macro MULTI_BUFFER_IO.
 *
 * Notes:
 *
 * 1) File and record length limitations
 *   The package was originally written for block-structured streaming
 *   devices like magnetic tape.  Such devices have hardware limits on
 *   the block size and maximum number of blocks on a single tape.  These
 *   limits are stored in environment variables TAPEIO_FBLOCKSIZE and
 *   TAPEIO_MAXFBLOCK for fixed-block media and TAPEIO_MAXVBLOCK for
 *   variable-block media.  There seems to be no limit placed on the
 *   number of blocks that can be stored on variable-block media.
 *   When the package was extended to disk files, a new i/o mode was
 *   introduced called transparent, which does unblocked transfers
 *   directly to/from the user's buffer.  For transparent i/o, the limits
 *   on record and file size are those of the underlying file system.
 *   For the sake of portability the following limits are assumed:
 *   			record_size < 1<<31 bytes
 *   			file_size < 1<<64 bytes
 *
 * Description : fixed block writing
 *
 * Written by Const Olchanski 11 Apr 1992, 6 Sept 1992,
 *                            26 Apr 1993
 *
 *
 */

#define MULTI_BUFFER_IO true
#define _FILE_OFFSET_BITS 64
#define MAX_JOURNAL_ENTRIES 0

#include <assert.h>
#include <tapeIO.h>
#include <ntypes.h>

#ifdef MULTI_BUFFER_IO
#include <pthread.h>
#endif

#ifdef USE_ALARM
#include <alarm.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>

#define MAX_TIOFILES 50

#define MODE_FREE  0         /* free slot        */
#define MODE_RD    1         /* fixed block read     */
#define MODE_WR    2         /* fixed block write    */
#define MODE_TRANSPARENT 3   /* transparent mode */
#define MODE_VRD   4         /* variable block read  */
#define MODE_VWR   5         /* variable block write */

#ifdef MULTI_BUFFER_IO

#define MULTI_BUFFER_SIZE 524288
#define MULTI_BUFFER_COUNT 3

typedef enum { empty, filled, flush, dummy } bufstat_t;

typedef struct buff_desc_
{
  bufstat_t status;      /* current status of this buffer */
  int size;		 /* size of this buffer in bytes */
  char *start;		 /* pointer to start of buffer area */
  int nbytes;		 /* byte count of current valid contents */
  int error;		 /* error code from most recent io op */
  off_t offset; 	 /* file offset to start of data (64 bits) */
  pthread_mutex_t mutex; /* pthreads mutex for thread synchronization */
  pthread_cond_t signal; /* pthreads cond variable for change of status */
  struct buff_desc_ *next; /* pointer to next buffer in list */
} buff_desc_t;

#endif

typedef struct 
{
  int fd;		 /* unix system file descriptor */
  FILE *fp;		 /* c library stream pointer (if any) */
  int ttype;		 /* type code for device */
  int mode;		 /* i/o mode for transfers */
  int debug;		 /* flag to request debug printout*/
  int blockSize;	 /* size of the internal buffer */
  int thisBlock;	 /* bytes in the current block */
  int inBuffer;		 /* unread bytes in the internal buffer */
  int rptr;		 /* read pointer (input descriptor only) */
  int maxnblocks;	 /* maximum size of file in blocks */
  off_t filemark; 	 /* file offset to valid data (64 bits) */
  uint64 bytecount;	 /* total i/o byte count (64 bits) */
  char *buffer;		 /* pointer to internal buffer */
#ifdef MULTI_BUFFER_IO
  buff_desc_t *bufdesc;	 /* pointer to buffer descriptor */
  pthread_t reader;      /* reader thread descriptor (if any) */
#endif
} tape_info_t;

static tape_info_t tape[MAX_TIOFILES];

static int nFd = 0;

static int getInfo(int fd);
static int setupTape(int fd,int mode);
static int read1(int fd,void*ptr,int sz);
static void *stream_reader(void *argp);
static int mbread(int p, char *buffer, int len);
static int mbflush(int p);

#ifdef MULTI_BUFFER_IO
static buff_desc_t *get_next_data_buffer(int p);
#endif

#if MAX_JOURNAL_ENTRIES

static char** journal = 0;
static int journal_next = 0;
static void add_journal(char *fmt, int datum);
static void print_journal();

static void add_journal(char *fmt, int datum)
{
   if (journal == 0) {
      int i;
      journal = (char**)malloc(MAX_JOURNAL_ENTRIES*sizeof(char*));
      for (i=0; i < MAX_JOURNAL_ENTRIES; i++) {
         journal[i] = malloc(101*sizeof(char));
         journal[i][0] = 0;
      }
   }
   sprintf(journal[journal_next],fmt,datum);
   journal_next = (journal_next + 1) % MAX_JOURNAL_ENTRIES;
}

static void print_journal()
{
   int i;
   for (i=0; i < MAX_JOURNAL_ENTRIES; i++) {
      printf("%s",journal[journal_next]);
      journal_next = (journal_next + 1) % MAX_JOURNAL_ENTRIES;
   }
}

#endif

static int getInfo(int fd)
{
  int i;
  for (i=0; i<nFd; i++)
    if (tape[i].mode != MODE_FREE)
      if (tape[i].fd == fd) return i;

/* new tape */

  return -1;
}

static int setupTape(int fd,int mode)
{
  int i,ttype;

  for (i=0; i<nFd; i++)
    if (tape[i].mode == MODE_FREE) break;

  if (i>=nFd)
    {
      if (nFd>=MAX_TIOFILES)
	{
	  fprintf(stderr,"tapeIO.setupTape: Error -  internal table overflow. Recompile with larger MAX_TIOFILES (current=%d)\n",MAX_TIOFILES);

#ifdef USE_ALARM
	  esend(ALARM_FATAL,"tapeIO.setupTape: internal table overflow","");
#endif
	  errno = 0;
	  return -1;
	}
      i = nFd;
      nFd++;
    }

  ttype = tape_devType(fd);
  if (ttype < 0)
    {
      fprintf(stderr,"tapeIO.setupTape: Notice - fd %d: Cannot determine device type, assuming STREAM\n",fd);
#ifdef USE_ALARM
      esend(ALARM_WARNING,"tapeIO_setupTape: Cannot determine device type, assuming STREAM","");
#endif
      ttype = TAPEIO_STREAM;
    }

  switch (ttype)
    {
    default:
      fprintf(stderr,"tapeIO.setupTape: Notice - fd %d: Bug: Unknown device type %d, assuming STREAM\n",fd,ttype);
#ifdef USE_ALARM
      esend(ALARM_WARNING,"tapeIO_setupTape: Unknown device type, assuming STREAM","");
#endif
      ttype = TAPEIO_STREAM;

    case TAPEIO_STREAM:  /* !!!!!!!!!!!!!!! careful there is no "break" above !!!!!!!!!!!!!!!! */ 
                         /* !!!!!!!!!!!!!!! "default" falls through here      !!!!!!!!!!!!!!!! */

      tape[i].mode = MODE_TRANSPARENT;
      tape[i].blockSize = 0;
      tape[i].maxnblocks = 0;
      break;

    case TAPEIO_FIXED:
      {
	char *sz = getenv("TAPEIO_FBLOCKSIZE");
	if (sz == NULL)
	  {
	    fprintf(stderr,"tapeIO.setupTape: Notice - the TAPEIO_FBLOCKSIZE environment variable is not set. Exit(1)...\n");
#ifdef USE_ALARM
	    esend(ALARM_FATAL,"tapeIO.setupTape: the TAPEIO_FBLOCKSIZE environment variable is not set","Exiting...");
#endif
	    exit(1);
	  }

	tape[i].blockSize = (int)strtoul(sz,NULL,0);

	if (tape[i].blockSize == 0)
	  {
	    fprintf(stderr,"tapeIO.setupTape: Error - the TAPEIO_FBLOCKSIZE envirnment variable (equal to (%s)) is zero or not a correct block size. Exit(1)...\n",sz);
#ifdef USE_ALARM
	    esend(ALARM_FATAL,"tapeIO.setupTape: the TAPEIO_FBLOCKSIZE environment variable is wrong","Exit(1)...");
#endif
	    exit(1);
	  }
	
	sz = getenv("TAPEIO_MAXFBLOCK");
	if (sz == NULL)
	  {
	    fprintf(stderr,"tapeIO.setupTape: Error - the TAPEIO_MAXFBLOCK environment variable is not set. Exit(1)...\n");
#ifdef USE_ALARM
	    esend(ALARM_FATAL,"tapeIO.setupTape: the TAPEIO_MAXFBLOCK environment variable is not set","Exit(1)...");
#endif
	    exit(1);
	  }

	tape[i].maxnblocks = (int)strtoul(sz,NULL,0) / tape[i].blockSize;

	if (tape[i].blockSize == 0)
	  {
	    fprintf(stderr,"tapeIO.setupTape: Error - the TAPEIO_MAXFBLOCK envirnment variable (equal to (%s)) is zero or not a correct block size. Exit(1)...\n",sz);
#ifdef USE_ALARM
	    esend(ALARM_FATAL,"tapeIO.setupTape: the TAPEIO_MAXFBLOCK environment variable is wrong","Exiting...");
#endif
	    exit(1);
	  }
	
	tape[i].mode = mode;

	fprintf(stderr,"tapeIO.setupTape: Error - fd %d: tape with block size %d bytes, max transfer size is %d blocks\n",
		fd,tape[i].blockSize,tape[i].maxnblocks);
	
      }
      break;

    case TAPEIO_VARIABLE:
      {
	char *sz = getenv("TAPEIO_MAXVBLOCK");
	if (sz == NULL)
	  {
	    fprintf(stderr,"tapeIO.setupTape: Notice - TAPEIO_MAXVBLOCK is not set. Using default: 102400 bytes\n");
#ifdef USE_ALARM
	    esend(ALARM_FATAL,"tapeIO_setupTape: the TAPEIO_MAXVBLOCK environment variable is not set","Using default...");
#endif
	    sz = "102400";
	  }

	tape[i].blockSize = (int)strtoul(sz,NULL,0);
	tape[i].maxnblocks = 1;

	if (mode==MODE_RD)
	  {
	    tape[i].blockSize *= 2;
	  }

	switch (mode)
	  {
	  case MODE_WR: tape[i].mode = MODE_VWR; break;
	  case MODE_RD: tape[i].mode = MODE_VRD; break;
	  }

	fprintf(stderr,"tapeIO.setupTape: Error - fd %d: tape with variable block size, max tape block size is %d bytes\n",
		fd,tape[i].blockSize);
      }

      break;
    }

  tape[i].fd = fd;
  tape[i].fp = 0;
  tape[i].ttype = ttype;

  if (getenv("TAPEIO_DEBUG"))
    tape[i].debug = (int)strtoul(getenv("TAPEIO_DEBUG"),NULL,0);
  else
    tape[i].debug = 0;

  tape[i].filemark = 0;
  tape[i].bytecount = 0;
  tape[i].inBuffer = 0;
  tape[i].rptr = 0;

  if (tape[i].blockSize > 0)
    {
      tape[i].buffer = malloc(tape[i].blockSize);

      if (!tape[i].buffer)
	{
	  fprintf(stderr,"tapeIO.setupTape: Error - fd %d: cannot allocate %d bytes for internal buffer, errno: %d (%s)\n",fd,tape[i].blockSize,errno,strerror(errno));
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"tapeIO_setupTape: cannot allocate memory for buffer","");
#endif
	  return -1;
	}
    }
  else tape[i].buffer = NULL;

  return i;
}

/* This external function was added to the tapeIO library to make it
 * possible to use the standard c io library functions (fread/fwrite)
 * instead of the unbuffered low-level system calls (read/write) that
 * were used in days of yore (eg. to write directly onto tape).  It
 * takes a stream pointer argument instead of the descriptor that is
 * used by all of the other tapeIO library functions.  For example,
 *
 *   fp = fopen("foobar.dat","w");
 *   tape_setStream(fp);
 *   ...
 *   tape_write(fileno(fp),buf,len);
 *   ...
 *   tape_flush(fileno(fp));
 *   fclose(fp);
 *
 * The call to tape_flush is proper, but not required, before fclose.
 * -rtj 11/2000
 */

// Replacement for standard fileno in clib
// that seems not to work with libpdcap
int fileno(FILE*fp)
{
  if (fp) {
     return fp->_fileno;
  }
  return -1;
}

int tape_setStream(FILE* fp)
{
  int p,fd;

  fd = fileno(fp);
  if (fd == -1) {
    fprintf(stderr,"tapeIO.tape_setStream: Error - invalid stream file fp.\n");
#ifdef USE_ALARM
    esend(ALARM_FATAL,"tapeIO.tape_setStream: invalid stream file fp","");
#endif
    errno = 0;
    return -1;
  }
  p = getInfo(fd);
  if (p<0)
    for (p=0; p<nFd; p++)
      if (tape[p].mode == MODE_FREE) break;

  if (p>=nFd)
    {
      if (nFd>=MAX_TIOFILES)
	{
	  fprintf(stderr,"tapeIO.tape_setStream: Error - internal table overflow.\n");
          fprintf(stderr,"Recompile with larger MAX_TIOFILES (current=%d)\n",MAX_TIOFILES);

#ifdef USE_ALARM
	  esend(ALARM_FATAL,"tapeIO.tape_setStream: internal table overflow","");
#endif
	  errno = 0;
	  return -1;
	}
      p = nFd;
      nFd++;
    }

  tape[p].ttype = TAPEIO_STREAM;
  tape[p].mode = MODE_TRANSPARENT;
  tape[p].blockSize = 0;
  tape[p].maxnblocks = 0;
  tape[p].fd = fd;
  tape[p].fp = fp;

  if (getenv("TAPEIO_DEBUG"))
    tape[p].debug = (int)strtoul(getenv("TAPEIO_DEBUG"),NULL,0);
  else
    tape[p].debug = 0;

  tape[p].filemark = 0;
  tape[p].bytecount = 0;
  tape[p].inBuffer = 0;
  tape[p].rptr = 0;
  tape[p].buffer = NULL;
#ifdef MULTI_BUFFER_IO
  tape[p].bufdesc = NULL;
#endif

  return 0;
}

FILE *tape_getStream(int fd)
{
  int p = getInfo(fd);
  if (p<0)
    return 0;
  else
    return tape[p].fp;
}

int tape_vwrite(int p)
{
  int err;

  if (tape[p].inBuffer & 1) /* for some wierd reason the blocksize should be even */
    {
      tape[p].buffer[tape[p].inBuffer] = 0;   /* pad with a zero */
      tape[p].inBuffer ++;
    }

  err = write(tape[p].fd,tape[p].buffer,tape[p].inBuffer);
  if (err<0)
    {
      fprintf(stderr,"tapeIO.tape_vwrite: Error - write() returned %d, errno: %d(%s)\n",err,errno,strerror(errno));
#ifdef USE_ALARM
      esend(ALARM_ERROR,"tapeIO.tape_vwrite: write() error","%s",strerror(errno));
#endif
    }

  tape[p].inBuffer = 0;
  
  return err;
}

int tape_vread(int p)
{
  int retry;

  tape[p].inBuffer = 0;
  tape[p].rptr = 0;

  do /* retry read() requests */
    {
      int rd;
      retry = 0;
      errno = 0;

      rd = read(tape[p].fd,tape[p].buffer,tape[p].blockSize);
      if (rd > 0)
	{
	  tape[p].thisBlock = rd;
	  tape[p].inBuffer = rd;

	  if (tape[p].debug)
	    fprintf(stderr,"tapeIO.tape_vread: Error - fd=%d, requested=%d, read=%d\n",tape[p].fd,tape[p].blockSize,rd);

	  retry = 0;
	}
      else
	{
	  switch (errno)
	    {
	    case 0:
	      switch (tape[p].ttype)
		{
		default:
		case TAPEIO_STREAM:
		  return TAPEIO_EOT;
		case TAPEIO_FIXED:
		case TAPEIO_VARIABLE:
		  return 0;            /* file mark   */
		}

	    case ENOSPC: return TAPEIO_EOT;   /* end of tape */
	    case ENOMEM:
	      /*
	       * tape block longer than our buffer,
	       * increase the buffer and retry does not work because the retry
	       * will give us the *next* block. This behaviour is not documented
	       * in the man pages (man read, mtio and tps), therefore the
	       * "move the tape back by one block and try to read again"
	       * would be unreliable as well.
	       *
	       * Easiest is to return an error and have somebody else worry about recovering.
	       *
	       * CO 30-Oct-1994.
	       */

	      fprintf(stderr,"tapeIO.tape_vread: Error - the next tape block cannot be read because it is longer than %d bytes, increase TAPEIO_MAXVBLOCK and try again.\n",tape[p].blockSize);
	      return -1;
	      
	    default:
	      fprintf(stderr,"tapeIO.tape_vread: Error - read() returned %d, errno: %d (%s)\n",rd,errno,strerror(errno));
#ifdef USE_ALARM
	      esend(ALARM_ERROR,"tapeIO.tape_vread: read() error","%s",strerror(errno));
#endif
	      return -1;
	    }
	}
    } while (retry);

/*printf("vread: %d, %d\n",tape[p].blockSize,tape[p].inBuffer);*/
  
  return tape[p].inBuffer;
}

int tape_flush(int fd)
{
  int err = 0;
  int p;
  int i;

  p = getInfo(fd);
  if (p<0) return 0;

  switch (tape[p].mode)
    {
    default:
    case MODE_RD:

#ifdef MULTI_BUFFER_IO
       if (tape[p].fp) {
          mbflush(p);
       }
#endif
       break;

    case MODE_WR:

      /* pad the last block with zeros */

      for (i=tape[p].inBuffer; i<tape[p].blockSize; i++)
	tape[p].buffer[i] = 0;

      /* write the last block onto tape */

      err = (int)write(tape[p].fd,tape[p].buffer,tape[p].blockSize);
      if (err<0)
	{
	  fprintf(stderr,"tapeIO.tape_flush: Error - write() returned %d, errno: %d (%s) (buffer flushed anyway)\n",err,errno,strerror(errno));
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"tapeIO.tape_flush: write() error"," %s (buffer flushed anyway)\n",strerror(errno));
#endif
	}
      break;

    case MODE_VRD:
      fprintf(stderr,"tapeIO.tape_flush: Notice - the input buffer was flushed: blocksize: %d, %d bytes in the current block, %d unread bytes were flushed\n",tape[p].blockSize,tape[p].thisBlock,tape[p].inBuffer);
      tape[p].inBuffer = 0;
      break;

    case MODE_VWR:
      err = tape_vwrite(p);
      break;
    }

  if (tape[p].buffer != 0) {
     free(tape[p].buffer);
     tape[p].buffer = NULL;
  }
  tape[p].mode = MODE_FREE;
  return err;
}

int tape_ff(int fd,off_t skipbytes)
{
  int p;
  int err = 0;
  off_t loc;

  p = getInfo(fd);
  if (p<0)
    {
      p = setupTape(fd,MODE_RD);
      if (p<0) return -1;
    }
  loc = (tape[p].bytecount > tape[p].filemark)?
         tape[p].bytecount : tape[p].filemark ;
  loc += skipbytes;

  switch (tape[p].mode)
    {
    default:
    case MODE_RD:

       if (tape[p].fp) {
#ifdef MULTI_BUFFER_IO
          buff_desc_t *this = tape[p].bufdesc;
          tape[p].filemark = loc;
#else
	  if ((err = fseeko(tape[p].fp,loc,SEEK_SET)) == 0) {
	     tape[p].bytecount += skipbytes;
	  }
#endif
	  return err;
       }
       else if ((err = lseek(fd,loc,SEEK_SET)) == loc) {
	  tape[p].bytecount += skipbytes;
	  return err = 0;
       }

    case MODE_WR:
    case MODE_VRD:
    case MODE_VWR:

        fprintf(stderr,"tapeIO.tape_ff: Error - ff operation allowed only in transparent mode\n");
        err = -1;
    }
  return err;
}

int tape_fwrite_fillBuf(int p,const void *ptr,int len)
{
  /* first fill the buffer */
	  
  int fill = tape[p].blockSize - tape[p].inBuffer;
	    
  if (len<fill) fill = len;
	    
  memcpy(&tape[p].buffer[tape[p].inBuffer],ptr,fill);
	    
  tape[p].inBuffer += fill;
  tape[p].bytecount += fill;

  /* now write out the buffer */
	  
  if (tape[p].inBuffer == tape[p].blockSize)
    {
      ssize_t err = write(tape[p].fd,tape[p].buffer,tape[p].inBuffer);
      if (err<0)
	{
	  fprintf(stderr,"tapeIO.tape_fwrite_fillBuf: Error - write1 error: write() returned %d, errno: %d (%s)\n",err,errno,strerror(errno));
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"tapeIO.fwrite_fillBuf: write() failed","%s\n",strerror(errno));
#endif
	  return -1;
	}
      tape[p].inBuffer = 0;
    }
  else
    {
      if (len==0) return fill;
      else 
	{
	  fprintf(stderr,"tapeIO.tape_fwrite_fillBuf: Error - Bug: This branch should never be executed. See source code (w1)\n");
	  /* explanation: buffer is not full, and yet len is not 0 */
	  errno = 0;
	  return -1;
	}
    }

  return fill;
}

int tape_write(int fd,const void*buffer,int len)
{
  int err = 0;
  int wr;
  int p;
  const char *ptr = buffer;

  if ((len & 3) != 0)
    {
      fprintf(stderr,"tapeIO.tape_write: Error - write call with invalid length: %d should be divisible by 4\n",len);
      return (-1);
    }

  p = getInfo(fd);
  if (p<0) 
    {
      p = setupTape(fd,MODE_WR);
      if (p<0) return -1;
    }

  switch (tape[p].mode)
    {
    default:
      fprintf(stderr,"tapeIO.tape_write: Error - write call with wrong file mode\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"tapeIO.tape_write: Write call with wrong file mode","");
#endif
      return -1;

    case MODE_TRANSPARENT:
/*
 * The standard c-library streams interface provides a more efficient
 * buffering facility than the primitive one implemented here.  Why
 * not let the stream do it for us, if it exists?
 */
      if (tape[p].fp) {
        wr = fwrite(buffer,1,len,tape[p].fp);
      }
      else {
        wr = (int)write(fd,buffer,len);
      }
      if (wr<0) {
	fprintf(stderr,"tapeIO.tape_write: Error - (direct) write() returned %d, errno: %d (%s)\n",wr,errno,strerror(errno));
#ifdef USE_ALARM
	esend(ALARM_ERROR,"tapeIO.tape_write: (direct) write() failed","%s",strerror(errno));
#endif
      }
      else
        tape[p].bytecount += wr;
      return wr;

    case MODE_VWR:
      if (tape[p].inBuffer + len >= tape[p].blockSize)
	err = tape_vwrite(p);

      if (tape[p].inBuffer + len >= tape[p].blockSize)   /* this event is too big to fit in one block */
	{
	  err = (int)write(tape[p].fd,buffer,len);
	  if (err<0)
	    {
	      fprintf(stderr,"tapeIO.tape_write: Error - write() returned %d, errno: %d (%s)\n",err,errno,strerror(errno));
#ifdef USE_ALARM
	      esend(ALARM_ERROR,"tapeIO.tape_write: write() error","%s",strerror(errno));
#endif
	      return -1;
	    }
	}
      else
	{
	  memcpy(&tape[p].buffer[tape[p].inBuffer],buffer,len);
	  tape[p].inBuffer += len;
	  err = len;
	}

      tape[p].bytecount += err;
      return err;

    case MODE_WR:
/*
      fprintf(stderr,"Got into MODE_WR...\n");
      abort();
*/
      wr = 0;

      if (tape[p].inBuffer>0) /* if buffer is empty, do not bother filling it */
	{
	  err = tape_fwrite_fillBuf(p,ptr,len);
	  if (err > 0)
	    {
	      len -= err;
	      ptr += err;
	      wr += err;
	    }
	}

      /* now write out data (if have enough of it) */
      
      {
	while (len>=tape[p].blockSize)
	  {
	    int count,ret;
	    int n = len/tape[p].blockSize;  /* number of complete blocks to write */
	    
	    if (n>tape[p].maxnblocks) n = tape[p].maxnblocks;
	    
	    count = n*tape[p].blockSize;
	    
	    ret = (int)write(tape[p].fd,ptr,count);
	    if (ret<0)
	      {
		fprintf(stderr,"tapeIO.tape_write: Error - (write2) write() returned %d, errno: %d (%s)\n",ret,errno,strerror(errno));
#ifdef USE_ALARM
		esend(ALARM_ERROR,"tapeIO.tape_write: write() failed","%s\n",strerror(errno));
#endif
		return -1;
	      }
	    ptr += ret;
	    len -= ret;
	    wr += ret;
	  }
      }
      
      /* put the rest into the buffer for next write or tape_flush() */
      
      {
	if (len>=tape[p].blockSize)
	  {
	    fprintf(stderr,"tapeIO.tape_write: Error - Bug: This branch should never be executed. See source code (w2)\n");
	    /* explanation: I was supposed to write all the full blocks I have */
	    errno = 0;
	    return -1;
	  }
	
	if (len>0)
	  {
	    memcpy(tape[p].buffer,ptr,len);
	    
	    tape[p].inBuffer = len;
	    ptr += len;
	    len -= len;
	    wr += len;
	  }
      }
      tape[p].bytecount += wr;
      return wr;

    }
}

int tape_read(int fd,void*buffer,int len)
{
  int ret;
  int rd;
  int p;
  char *ptr = buffer;

  p = getInfo(fd);
  if (p<0) 
    {
      p = setupTape(fd,MODE_RD);
      if (p<0) return -1;
    }
/*
fprintf(stderr,"tapeIO.tape_read: %d, 0x%x, %d\n",fd,buffer,len);
*/
  switch (tape[p].mode)
    {
    default:
      fprintf(stderr,"tapeIO.tape_read: Error - read call with wrong file mode\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"tapeIO.tape_read: read call with wrong file mode","");
#endif
      return -1;

    case MODE_TRANSPARENT:

      if (tape[p].fp) {
#ifdef MULTI_BUFFER_IO
        rd = mbread(p,(char*)buffer,len);
#else
        rd = fread(buffer,1,len,tape[p].fp);
#endif
      }
      else {
        rd = read1(fd,buffer,len);
      }
      if (rd<0) {
	fprintf(stderr,"tapeIO.tape_read: Error - (direct) read() returned %d, errno: %d (%s)\n",rd,errno,strerror(errno));
#ifdef USE_ALARM
	esend(ALARM_ERROR,"tapeIO.tape_read: (direct) read() failed","%s",strerror(errno));
#endif
      }
      else if (rd==0)
	return TAPEIO_EOT;
      else
        tape[p].bytecount += rd;

      return rd;

    case MODE_VRD:

      if (tape[p].inBuffer < 2) /* 0 or 1 - 0 means buffer exhausted, 1 - the buffer was padded to even bytes */
	{
	  ret = tape_vread(p);
	  if (ret <= 0) return ret;
	}

      if (len > tape[p].inBuffer)
	{
	  fprintf(stderr,"tapeIO.tape_read: Error - attempt to read %d bytes when buffer contains only %d bytes, blocksize is %d, this block has %d bytes\n",len,tape[p].inBuffer,tape[p].blockSize,tape[p].thisBlock);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"tapeIO.tape_read: attempt to read more bytes then in buffer","%d vs %d",len,tape[p].inBuffer);
#endif
	  return -1;
	}

      memcpy(buffer,&tape[p].buffer[tape[p].rptr],len);
      tape[p].rptr += len;
      tape[p].inBuffer -= len;
      tape[p].bytecount += len;
      return len;

    case MODE_RD:

      rd = 0;

      if (tape[p].inBuffer>0) /* if buffer not empty, read from it */
	{
	  
	  /* first read from buffer */
	  
	  {
	    int fill = tape[p].inBuffer;
	    
	    if (len<fill) fill = len;
	    
	    memcpy(ptr,&tape[p].buffer[tape[p].rptr],fill);
	    
	    len -= fill;
	    ptr += fill;
	    tape[p].rptr += fill;
	    tape[p].inBuffer -= fill;
	    rd += fill;
	  }
	  
	}

      /* sanity check */
      
      if (len>0)
	if (tape[p].inBuffer!=0)
	  {
	    fprintf(stderr,"tapeIO.tape_read: Error - Bug: This branch should never be executed. See source code (r1)\n");
	    /* explanation: if len is non zero then all of the buffer should be read out */
	    errno = 0;
	    return -1;
	  }

      /* now read data from input stream */
      
      {
	if (len>=tape[p].blockSize)
	  {
	    int n = len/tape[p].blockSize;
	    int count = n*tape[p].blockSize;
	    
	    ret = (int)read(tape[p].fd,ptr,count);
	    if (ret<0)
	      {
		fprintf(stderr,"tapeIO.tape_read: Error - read() (@data) returned %d, errno: %d (%s)",ret,errno,strerror(errno));
#ifdef USE_ALARM
		esend(ALARM_ERROR,"tapeIO.tape_read: read() failed","%s",strerror(errno));
#endif
		return -1;
	      }
	    
	    if (ret==0) return rd;  /* eof detected */
	    
	    if (ret!=count) {
	      fprintf(stderr,"tapeIO.tape_read: Error - Short read() to @data, read(%d) returned %d, errno: %d (%s)",count,ret,errno,strerror(errno));
#ifdef USE_ALARM
	      esend(ALARM_ERROR,"tapeIO.tape_read: Short read()","%s",strerror(errno));
#endif
	    }

	    ptr += ret;
	    len -= ret;
	    rd += ret;
	  }
      }

      /* if have more to read but it is less than a full block, read it into buffer */

      {
	if (len>=tape[p].blockSize)
	  {
	    fprintf(stderr,"tapeIO.tape_read: Error - Bug: This branch should never be executed. See source code (r2)\n");
	    /* explanation: I was supposed to read all the full blocks I have */
	    errno = 0;
	    return -1;
	  }

	if (len>0)
	  {
	    ret = (int)read(tape[p].fd,tape[p].buffer,tape[p].blockSize);
	    if (ret<0)
	      {
		fprintf(stderr,"tapeIO.tape_read: Error - read(%d) (@buffer) returned %d, errno: %d (%s)\n",tape[p].blockSize,ret,errno,strerror(errno));
#ifdef USE_ALARM
		esend(ALARM_ERROR,"tapeIO.tape_read: read() failed","%s",strerror(errno));
#endif
		return -1;
	      }
	    
	    tape[p].inBuffer = ret;
	    tape[p].rptr = 0;

	    if (ret==0)    /* end of file */
	      {
                tape[p].bytecount += rd;
		return rd;
	      }
	    
	    if (len>tape[p].inBuffer)
	      {
		fprintf(stderr,"tapeIO.tape_read: Error - Bug: This branch should never be executed. See source code (r3)\n");
		/* explanation: all the full blocks are expected to be read before */
		errno = 0;
		return -1;
	      }

	    memcpy(ptr,&tape[p].buffer[tape[p].rptr],len);

	    tape[p].inBuffer -= len;
	    tape[p].rptr += len;
	    ptr += len;
	    rd += len;
	    len -= len;
	    
	    if (len!=0)
	      {
		fprintf(stderr,"tapeIO.tape_read: Error - short read()? to buffer\n");
#ifdef USE_ALARM
		esend(ALARM_ERROR,"tapeIO.tape_read: Short read() to buffer","%s",strerror(errno));
#endif
	      }
	  }
      }
      
      tape[p].bytecount += rd;
      return rd;

    } /* end switch */
}

static int read1(int fd,void*ptr,int sz)
{
  char*p=ptr;
  int rd,ret;
  rd = 0;
  while (sz>0)
    {
      ret = (int)read(fd,p,sz);
      if (ret<0) return ret;
      if (ret==0) return rd;
      rd += ret;
      sz -= ret;
      p += ret;
    }
  return rd;
}

#ifdef MULTI_BUFFER_IO

static void *stream_reader(void *argp)
{
   static int one=1;
   int p = *(int *)argp;
   FILE *fp = tape[p].fp;
   buff_desc_t *this = tape[p].bufdesc;
   off_t stream_offset = ftello(fp);
   off_t fastforward_offset;
   free(argp);

   assert (this != 0);     /* the buffer list should already be there */

   this = this->next;
   pthread_mutex_lock(&this->mutex);
   pthread_cond_signal(&this->signal);

   while (fp) {
      this->offset = stream_offset;
      fastforward_offset = tape[p].filemark;
      if (fastforward_offset > stream_offset) {
         if (fseeko(fp,fastforward_offset,SEEK_SET)) {
            this->error = ferror(fp);
         }
         else {
            this->error = 0;
            stream_offset = fastforward_offset;
	    this->nbytes = stream_offset - this->offset;
         }
         this->status = dummy;
      }

      if (this->status == empty) {
         stream_offset += this->nbytes = fread(this->start,1,this->size,fp);
         if (this->nbytes != this->size) {
            this->error = ferror(fp);
         }
         else {
            this->error = 0;
         }
         this->status = filled;
      }

      pthread_mutex_lock(&this->next->mutex);
      pthread_mutex_unlock(&this->mutex);
      this = this->next;
      if (this->status != empty) {
         pthread_mutex_unlock(&this->mutex);
         break;
      }
   }
   return &one;
}

static buff_desc_t *get_next_data_buffer(int p)
{
   int *argp;
   FILE *fp = tape[p].fp;
   buff_desc_t *this = tape[p].bufdesc;
   if (this == 0) {
      buff_desc_t **link = &tape[p].bufdesc;
      int i;
      for (i=0; i<MULTI_BUFFER_COUNT; i++) {
         this = malloc(sizeof(buff_desc_t));
         this->status = empty;
         this->size = MULTI_BUFFER_SIZE;
         this->start = malloc(MULTI_BUFFER_SIZE);
         pthread_mutex_init(&this->mutex, NULL);
         pthread_cond_init(&this->signal, NULL);
         this->nbytes = 0;
         this->error = 0;
	 *link = this;
	 link = &this->next;
      }
      *link = this = tape[p].bufdesc;
      pthread_mutex_lock(&this->mutex);
      pthread_mutex_lock(&this->next->mutex);
      argp = malloc(sizeof(int));
      *argp = p;
      pthread_create(&tape[p].reader, NULL,
                     &stream_reader, (void *)argp);
      pthread_cond_wait(&this->next->signal,&this->next->mutex);
      pthread_mutex_unlock(&this->next->mutex);
   }

   pthread_mutex_lock(&this->next->mutex);
   this->status = empty;
   pthread_mutex_unlock(&this->mutex);
   this = this->next;

   assert (this->status != empty);
 
#if MAX_JOURNAL_ENTRIES
   add_journal("get_next_data_buffer:\n",0);
   add_journal("    status: %d\n",this->status);
   add_journal("    bytes: %d\n",this->nbytes);
   add_journal("    error: %d\n",this->error);
   add_journal("    offset: %d\n",this->offset);
#endif

   return tape[p].bufdesc = this;
}

static int mbread(int p, char *buffer, int len)
{
   int count;
   char *dst = buffer;
   buff_desc_t *this = tape[p].bufdesc;
#if MAX_JOURNAL_ENTRIES
   add_journal("mbread entry:\n",0);
   add_journal("   rptr = %d\n",tape[p].rptr);
   add_journal("   inBuffer = %d\n",tape[p].inBuffer);
   add_journal("   bytecount = %d\n",tape[p].bytecount);
   add_journal("   filemark = %d\n",tape[p].filemark);
#endif
   if (this == 0 || tape[p].inBuffer <= 0) {
      this = get_next_data_buffer(p);

      assert(this != NULL);   /* get_next_data_buffer should never
                                 return empty-handed!              */
      if (this->error) {
         fprintf(stderr,"tapeIO.mbread: Error - read returned %d, errno: %d (%s)\n",
                 this->size, this->error, strerror(this->error));
         return -1;
      }
      tape[p].inBuffer  = this->nbytes;
      tape[p].blockSize = this->size;
      tape[p].rptr      = 0;
   }
#if MAX_JOURNAL_ENTRIES
   add_journal("mbread checkpoint 1:\n",0);
   add_journal("   rptr = %d\n",tape[p].rptr);
   add_journal("   inBuffer = %d\n",tape[p].inBuffer);
   add_journal("   bytecount = %d\n",tape[p].bytecount);
   add_journal("   filemark = %d\n",tape[p].filemark);
#endif
   if (tape[p].filemark > tape[p].bytecount) {
      off_t skipbytes = tape[p].filemark - tape[p].bytecount;
      skipbytes -= tape[p].inBuffer;
      while (skipbytes >= 0) {
        this = get_next_data_buffer(p);
        if (this->error) {
          fprintf(stderr,"tapeIO.mbread: Error - read returned %d, errno: %d (%s)\n",
                  this->size, this->error, strerror(this->error));
          return -1;
        }
        else if (this->nbytes == 0) {
          skipbytes = 0;
	  break;
        }
        skipbytes -= this->nbytes;
      }
      tape[p].blockSize = this->size;
      tape[p].inBuffer = -skipbytes;
      tape[p].rptr = this->nbytes + skipbytes;
      tape[p].bytecount = tape[p].filemark;
      tape[p].filemark = 0;
   }
#if MAX_JOURNAL_ENRIES
   add_journal("mbread checkpoint 2:\n",0);
   add_journal("   rptr = %d\n",tape[p].rptr);
   add_journal("   inBuffer = %d\n",tape[p].inBuffer);
   add_journal("   bytecount = %d\n",tape[p].bytecount);
   add_journal("   filemark = %d\n",tape[p].filemark);
#endif
   if (tape[p].inBuffer <= 0) {
      return 0;
   }
   else if (this->status != filled) {
      fprintf(stderr,"tapeIO.mbread: Error - bad buffer found in pool.\n");
#if MAX_JOURNAL_ENTRIES
      print_journal();
#endif
      return -1;
   }
   count = (len < tape[p].inBuffer) ? len : tape[p].inBuffer;
   memcpy(dst,&tape[p].bufdesc->start[tape[p].rptr],count);
   tape[p].inBuffer -= count;
   tape[p].rptr += count;
   dst += count;
   len -= count;
   if (len > 0) {
      int rem = mbread(p,dst,len);
      return (rem >= 0)? count+rem : rem;
   }
   return count;
}

static int mbflush(int p)
{
   int i;
   void *retval;
   buff_desc_t *this = tape[p].bufdesc;
   if (this == 0) {
      return 0;
   }
   this->status = flush;
   pthread_mutex_unlock(&this->mutex);
   pthread_join(tape[p].reader, &retval);
   for (i=0; i<MULTI_BUFFER_COUNT; i++) {
      buff_desc_t *next = this->next;
      pthread_mutex_destroy(&this->mutex);
      pthread_cond_destroy(&this->signal);
      free(this->start);
      free(this);
      this = next;
   }
   tape[p].bufdesc = 0;
   return 0;
}

#endif


/* end file */
