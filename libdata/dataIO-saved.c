/*
 * dataIO.c
 *
 * Modification history
 *
 * July 30, 2004 - Richard Jones
 * - corrected a bug where a sieve file written with this library could
 *   not be read by it.  The situation occurred when a sieve contained a
 *   reference to a file but wanted to skip 100% of its contents.  The
 *   sieve semantics for this case is a single line containing a single 0
 *   after the input="" line.  On input this produced a "sieve overrun"
 *   error, which means that the sieve byte counts are chopping off the
 *   end of an event record.  The solution was to rearrange the sieve-open
 *   logic at the start data_read so that read pointers are always advanced
 *   to the start of the next valid data before the read passes down to 
 *   the file descriptors referenced by the sieve.
 *
 * February 22, 2003 - Richard Jones
 * - added new indirection mechanism for data files called "sieves".
 *   Sieve files function like symbolic links with additional filter
 *   functionality.  They do not contain the actual data but hold a
 *   symbolic link to the source file together with a binary mask that
 *   can be used to select arbitrary portions of the source file in the
 *   replica.  The order of the input data is not changed but arbitrary
 *   pieces of the file can be turned on or off.  From the user's point
 *   of view, reading from a sieve file is exactly the same as reading
 *   from a copy of the source file with the skipped portions compressed
 *   out.  Like symbolic links, sieve files can point to other sieve files,
 *   with cumulative effects.  To access a sieve file a new function
 *   data_setSieve must be invoked after the input and output files have
 *   been opened.
 *
 * February 17, 2003 - Richard Jones
 * - introduced a new function data_ff to provide a mechanism to skip
 *   quickly over unwanted regions of an input file.  It is implemented
 *   with immediate return; any wait is incurred upon the following read.
 *
 *
 * Notes:
 * 1. How to create a sieve file.
 *   A data sieve is an efficient means for saving a subset of data from
 *   a source file in an output file without creating a new copy of the
 *   data.  Once an output file is designated as a sieve to an input file,
 *   calls to data_write on that unit no longer store the data themselves
 *   in the output file but instead write a record of the regions of the
 *   input file that were selected.  The selected region is that segment
 *   of the input file that was last returned by a call to data_read.
 *   The way sieve files are created is shown in the following example.
 *
 *   FILE *fp,*fpin;
 *   fp = fopen("foobar.dat","w");
 *                                       <--- start of repeat block
 *   *fpin = fopen("foo.dat","r");
 *   int fd = fileno(fp);
 *   int fdin = fileno(fpin);
 *   data_setSieve(fp,fdin,"foo.dat");
 *   ...
 *   while (data_read(fdin,buf,len)) {
 *      if (...) {
 *         data_write(fd,buf,len);
 *      }
 *   }
 *   ...
 *   data_flush(fd);
 *   fclose(fpin);
 *                                       <--- end of repeat block
 *   fclose(fp);
 *
 *   Repetitions of the above repeat block will produce a sieve file that
 *   combines any number of source data files into one stream.  The call
 *   to data_flush is required when the input file is closed to allow the
 *   sieve record to be properly terminated.  Sieve data may not be mixed
 *   with ordinary data in the same output file; it is possible to write
 *   such a file but trying to read from it will produce an error.  Note
 *   that calling data_setSieve with fdin=0 sets up fp for reading instead
 *   of writing; sieving is not allowed on standard input.
 *
 * 2. How to read from a sieve file.
 *   Ideally reading from a sieve file should be indistinguishable from
 *   reading from an ordinary data file from the user's point of view.
 *   In fact, one extra step is required to declare the input as a sieve
 *   file BEFORE attempting to read from it, as follows.
 *
 *   FILE *fp = fopen("foobar.dat","r");
 *   int fd = fileno(fp);
 *   data_setSieve(fp,0,0);
 *   ...
 *   while (data_read(fdin,buf,len)) {
 *      ...
 *   }
 *   ...
 *   fclose(fp);
 *
 *   Invoking data_setSieve before reading on fp is necessary because
 *   the tapeIO interface which is ordinarily used to read and write
 *   buffered stream data is inefficient for line-oriented text.  The
 *   setSieve call enables data_read to bypass the tapeIO layer and
 *   directly access the sieve data through the streams interface.
 *   If the file is opened with evt_open() then the data_setSieve call
 *   is made automatically on sieve files so the sieve is transparent.
 *
 * 3. The format of a sieve data stream
 *   Sieve data files are written entirely in ascii.  This makes them
 *   easy to edit and transport, and does not cost much in terms of
 *   data volume.  It also makes them insensitive to byte-order issues.
 *   In the following diagram the vertical bar delimits each line.
 *   |DATA SIEVE version 1.0
 *   |input="/rooted/path/to/file.itape"
 *   |r0,s0
 *   |r1,s1
 *   |r2,s2
 *   |...
 *   |rn
 *   The numeric strings r0,r1... and s0,s1... are the byte-lengths of
 *   segments of the input stream to read (r0,r1...) or skip (s0,s1...).
 *   All values ri and si are non-negative integers and are listed in
 *   pairs except for the last r value, which appears alone on a line.
 *   Any number of such DATA SIEVE sections may appear in a sieve file.
 *
 * 4. Restrictions on use of sieve files
 *   Sieve files are similar to symbolic links in that if the location
 *   or name of the original data file is changed they become invalid. 
 *   The sieve file itself is not endian-sensitive and can be written
 *   and read on machines of opposite byte-order.
 */

/* Depending on the sparseness of the sieve,
 * the following may help or hurt performance;
 * use STREAM_THROUGH_SIEVES with caution.
 */

#define STREAM_THROUGH_SIEVES true

#define _FILE_OFFSET_BITS 64
#define SIEVE_VERSION "1.0"

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#include <disData.h>
#include <dataIO.h>
#include <tapeIO.h>
#include <crc.h>
#include <ntypes.h>

#ifdef USE_ALARM
#include <alarm.h>
#endif

#define WSIZE 4      /* word size on tape (length word, CRC32 word, etc... */

typedef struct {
  FILE *fp;		 /* stream pointer to this descriptor */
  int sieve;		 /* flag indicating sieve (1) or not (0) */
  int link;		 /* index of linked descriptor (sieve only) */
  int writing;		 /* flag indicating writing (1) or reading (0) */
  int wrongendian;       /* flag indicating data is stored byte-swapped */
  uint64 lastcount; 	 /* byte count of last op (64 bits) */
  uint64 leftcount;  	 /* leftover from last op (64 bits) */
  uint64 bytecount;	 /* total i/o byte count (64 bits) */
} file_info_t;

static file_info_t *file_info[FOPEN_MAX];

static file_info_t *getInfo(int fd)
{
  static int initialized = 0;
  if (! initialized) {
    int i;
    for (i=0;i<FOPEN_MAX;i++) {
      file_info[i] = 0;
    }
    initialized = 1;
  }
  else if (fd >= FOPEN_MAX) {
     fprintf(stderr,"dataIO.getInfo: Error - too many files open at once.\n");
     exit(1);
  }
  return file_info[fd];
}

int data_setSieve(FILE *fp, int fdin, char* filename)
{
  int fd;
  file_info_t *infd, *infdin;

  fd = fileno(fp);
  if (fd == -1) {
    fprintf(stderr,"dataIO.data_setSieve: Error - input file is not open.\n");
    errno = 0;
    return -1;
  }
  if ((infd=getInfo(fd)) == 0) {
    file_info[fd] = infd = malloc(sizeof(file_info_t));
  }
  infd->fp = fp;
  infd->sieve = 0;
  if (fdin) {
    infd->link = fdin;
    infd->writing = 1;
  }
  else {
    infd->link = 0;
    infd->writing = 0;
  }
  infd->lastcount = 0;
  infd->leftcount = 0;
  infd->bytecount = 0;

  if (fdin) {
    if ((infdin=getInfo(fdin)) == 0) {
      file_info[fdin] = infdin = malloc(sizeof(file_info_t));
      infdin->sieve = 0;
      infdin->link = 0;
    }
    infdin->writing = 0;
    infdin->lastcount = 0;
    infdin->leftcount = 0;
    infdin->bytecount = 0;   /* assumed to be at start of stream */
  }

  /* make sure that fdin matches filename */

  if (fdin || filename) {
    struct stat fbuf,sbuf;
    if (fstat(fdin,&fbuf)==0 && stat(filename,&sbuf)==0) {
      if (strncmp((char *)&fbuf,(char *)&sbuf,sizeof(struct stat))) {
        fdin = 0;
      }
    }
    else {
      fdin = 0;
    }
    if (fdin == 0) {
      fprintf(stderr,"dataIO.data_setSieve: Error - input descriptor and filename do not match.\n");
      errno = 0;
      return -1;
    }
  }

  if (infd->writing) {
    fprintf(fp,"DATA SIEVE version %s\n",SIEVE_VERSION);
    fprintf(fp,"input=\"%s\"\n",filename);
  }
  infd->sieve = 1;
}

int data_write(int fd,const void* buffer_)
{
  int ret,size;
  void* buffer = (void*)buffer_;

  file_info_t *infd = getInfo(fd);
  if (infd == 0) {
    file_info[fd] = infd = malloc(sizeof(file_info_t));
    infd->fp = tape_getStream(fd);
    infd->sieve = 0;
    infd->link = 0;
    infd->writing = 1;
    infd->lastcount = 0;
    infd->leftcount = 0;
    infd->bytecount = 0;
  }
  else if (infd->sieve && infd->writing) {
    int wr = 0;
    int fdin = infd->link;
    if (fdin) {
      file_info_t *infdin = getInfo(fdin);
      const uint64 maxi = 1LL<<60;
      uint64 step,drop,load;
      step = infdin->bytecount - infd->bytecount;
      load = infd->bytecount + infd->leftcount;
      drop = step - infdin->lastcount;
      while (drop > 0 || load > maxi) {
        uint64 r = infd->leftcount;
        uint64 s = (drop < maxi)? drop : maxi;
        wr = fprintf(infd->fp,"%lld,%lld\n",r,s);
        infd->bytecount += s;
        infd->leftcount = 0;
        load -= r;
        drop -= s;
      }
      if (wr >= 0) {
        wr = infdin->lastcount;
	infd->leftcount += wr;
	infd->bytecount += wr;
	infd->lastcount = wr;
      }
    }
    else {
      fprintf(stderr,"dataIO.data_write: Error - write to sieve after data_flush.\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_write: write to sieve after data flush.\n");
#endif
      return -1;
    }
    return DATAIO_OK;
  }

  size = *(int*)buffer;
  if (size % WSIZE)
    {                                /* force size to be multiple of 4 (WSIZE) */
      size = ((size/WSIZE)+1)*WSIZE;
/*
fprintf(stderr,"Change size: %d -> %d\n", *(uint32*)buffer,size);
*/
      *(uint32*)buffer = size;
    }

  data_addCRC(buffer);

  ret = tape_write(fd,buffer,size + WSIZE);
  if (ret<0)
    {
      fprintf(stderr,"dataIO.data_write: Error - data write failure: tape_write() failed\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.write: data write failure","tape_write() failed\n");
#endif
      return DATAIO_ERROR;
    }
  else if (infd)
    {
      infd->bytecount += ret;
      infd->lastcount = ret;
    }

  return DATAIO_OK;
}

static uint32 data_swapEndianInt(uint32 w)
{
  return ((w&0xff)<<24)|((w&0xff00)<<8)|((w&0xff0000)>>8)|((w&0xff000000)>>24);
}
  
int data_read(int fd,void* buffer,int bufsize)
{
  int rd = 0;
  int ret, rec;
  char *buf = buffer;
  uint32 lsize; /* size in local-endian format */
  uint32 tsize; /* size on tape (maybe endian-swapped) */
  uint32 reclen;

  file_info_t *infd = getInfo(fd);
  if (infd == 0) {
    file_info[fd] = infd = malloc(sizeof(file_info_t));
    infd->fp = tape_getStream(fd);
    infd->sieve = 0;
    infd->link = 0;
    infd->writing = 0;
    infd->lastcount = 0;
    infd->leftcount = 0;
    infd->bytecount = 0;
    infd->wrongendian = 0;
    if (buffer == 0) {
      return DATAIO_OK;
    }
  }
  else if (infd->sieve) {
    if (infd->link == 0) {
      char filename[999];
      char version[99];
      FILE *fpin;
      char str[80];
      int hdr = sscanf(fgets(str,80,infd->fp),
                       "DATA SIEVE version %s\n",version);
      if (hdr == 0) {
        fprintf(stderr,"dataIO.data_read: Error - sieve header format error\n");
#ifdef USE_ALARM
        esend(ALARM_ERROR,"dataIO.data_read: sieve header format error\n");
#endif
        return DATAIO_ERROR;
      }
      else if (hdr == EOF) {
        return DATAIO_EOF;
      }
      else if ((sscanf(fgets(str,80,infd->fp),"input=\"%s\"\n",filename) == 1)
               && (fpin = fopen(strtok(filename,"\""),"r"))) {
        int is_nested_sieve = sscanf(fgets(str,80,fpin),
                                     "DATA SIEVE version %s\n",version);
        rewind(fpin);
        if (is_nested_sieve) { 
          data_setSieve(fpin,0,0);
	}
#if defined STREAM_THROUGH_SIEVES
	else {
          tape_setStream(fpin);
	}
#endif
	if (sscanf(fgets(str,80,infd->fp),"%lld",&infd->leftcount) != 1) {
          fprintf(stderr,"dataIO.data_read: Error - bad sieve data.\n");
#ifdef USE_ALARM
          esend(ALARM_ERROR,"dataIO.data_read: bad sieve data.\n");
#endif
	  return DATAIO_ERROR;
	}
        infd->link = fileno(fpin);
	if (infd->link == -1) {
          fprintf(stderr,"dataIO.data_read: Error - fpin is not valid.\n");
#ifdef USE_ALARM
          esend(ALARM_ERROR,"dataIO.data_read: fpin is not valid.\n");
#endif
	  return DATAIO_ERROR;
	}
        data_read(infd->link,0,0);
        while (infd->leftcount == 0) {
          if (data_read(fd,0,0) != DATAIO_OK) {
            break;
          }
        }
      }
      else {
        fprintf(stderr,"dataIO.data_read: Error - cannot open sieve input\n");
#ifdef USE_ALARM
        esend(ALARM_ERROR,"dataIO.data_read: cannot open sieve input\n");
#endif
	perror(filename);
        return DATAIO_ERROR;
      }
    }

    if (buffer) {
      if (infd->leftcount == 0) {
        return DATAIO_EOF;
      }

      ret = data_read(infd->link,buffer,bufsize);
    
      if (ret == DATAIO_OK) {
        int count = *(int*)buffer + WSIZE;
        if (count > infd->leftcount) {
          fprintf(stderr,"dataIO.data_read: Error - sieve data overrun.\n");
#ifdef USE_ALARM
          esend(ALARM_ERROR,"dataIO.data_read: sieve data overrun.\n");
#endif
          ret = DATAIO_ERROR;
	}
	else {
          infd->bytecount += count;
          infd->lastcount = count;
          infd->leftcount -= count;
	}
      }
    }
    else {
      ret = DATAIO_OK;
    }

    while (infd->leftcount == 0) {
      off_t skipbytes;
      char str[80];
      int cont = sscanf(fgets(str,80,infd->fp),",%lld\n",&skipbytes);
      if (cont == 0) {
        sscanf(fgets(str,80,infd->fp),"\n");
        data_flush(fd);
	break;
      }
      else if (cont != 1) {
        fprintf(stderr,"dataIO.data_read: Error - bad sieve data.\n");
#ifdef USE_ALARM
        esend(ALARM_ERROR,"dataIO.data_read: bad sieve data.\n");
#endif
        return DATAIO_ERROR;
      }
      if (skipbytes > 0) {
        data_ff(infd->link,skipbytes);
      }
      if (sscanf(fgets(str,80,infd->fp),"%lld",&infd->leftcount) != 1) {
        fprintf(stderr,"dataIO.data_read: Error - bad sieve data.\n");
#ifdef USE_ALARM
        esend(ALARM_ERROR,"dataIO.data_read: bad sieve data.\n");
#endif
        return DATAIO_ERROR;
      }
    }
    return ret;
  }

  ret = tape_read(fd,&tsize,WSIZE);

  if (ret==TAPEIO_EOT)
    return DATAIO_EOT;

  if (ret<0)
    {
      fprintf(stderr,"dataIO.data_read: Error - length word read failure: tape_read() failed\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_read: length word read failure","tape_read() failed\n");
#endif
      return DATAIO_ERROR;
    }
  else
    {
      rd += WSIZE;
    }
/*
fprintf(stderr,"dataIO.read: %d, 0x%x, %d, size: %d\n",fd,buffer,bufsize,size);
*/
  if (ret==0) return DATAIO_EOF;

  lsize = tsize;
  if (infd->wrongendian)
    {
      lsize = data_swapEndianInt(tsize);
    }

  /* fprintf(stderr,"Lsize: 0x%x 0x%x 0x%x\n",tsize,lsize,data_swapEndianInt(tsize)); */

  if (lsize>bufsize)
    {
      /* bogus size word: maybe we are the wrong-endian? */

      lsize = data_swapEndianInt(lsize);

      if (lsize < bufsize)
	{
	  /* now it's better. Assume swapped-endian and go ahead! */
	  
	  infd->wrongendian = !(infd->wrongendian);

	  fprintf(stderr,"dataIO.data_read: Notice - itape is wrong-endian, converting on-the-fly.\n");
	}
      else
	{
	  fprintf(stderr,"dataIO.data_read: Error -  record length (%d) is bigger than buffer (%d)\n",lsize,bufsize);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"dataIO.read: record length is bigger than buffer","(%d) vs. (%d)\n",lsize,bufsize);
#endif
	  return DATAIO_ERROR;
	}
    }

  if (lsize < WSIZE)
    {
      fprintf(stderr,"dataIO.data_read: Error - record length %d is less than %d\n",lsize,WSIZE);
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_read: record length too small","%d vs. %d\n",lsize,WSIZE);
#endif
      return DATAIO_ERROR;
    }

  *((uint32*)buf) = tsize;  /* use the tape size- it will be endian-swapped later */

  reclen = lsize - WSIZE + WSIZE; /* -WSIZE is the length word,
				     +WSIZE is the CRC32 word */
  rec = tape_read(fd,buf + WSIZE, reclen);
  if (rec<0)
    {
      fprintf(stderr,"dataIO.data_read: Error - data read failure: tape_read() failed\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_read: data read failure","tape_read() failed\n");
#endif
      return DATAIO_ERROR;
    }
  else if (rec!=reclen)
    {
      fprintf(stderr,"dataIO.data_read: Error - short read: %d bytes expected, %d bytes read\n",reclen,rec);
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_read: short read","event truncated\n");
#endif
      return DATAIO_ERROR;
    }
  else
    {
      rd += rec;
    }

  if (data_checkCRC1(buffer,lsize))  /* CRC32 is endian-neutral */
    {
      fprintf(stderr,"dataIO.data_read: Warning - data failed the CRC check\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_read: data failed the CRC check","");
#endif
      return DATAIO_BADCRC;
    }

  if (infd->wrongendian)
    endian_convertItape(buffer,lsize);

  infd->bytecount += rd;
  infd->lastcount = rd;
  return DATAIO_OK;
}

#define MAXBUF (512*1024)  /* 512 Kbytes */

int data_read_alloc(int fd,void **buffer)
{
  int ret,rd;
  char *buf = NULL;
  uint32 tsize;
  uint32 lsize;

  if (buffer==NULL)
    {
      fprintf(stderr,"dataIO.data_read_alloc: Error - NULL pointer\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_read_alloc: NULL buffer pointer","");
#endif
      return DATAIO_ERROR;
    }

  buf = malloc(lsize + WSIZE);

  if (!buf)
    {
      fprintf(stderr,"dataIO.data_read_alloc: Error - cannot allocate %d bytes for data buffer, errno: %d (%s)\n",lsize + WSIZE,errno,strerror(errno));
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_read_alloc: cannot allocate memory","%d bytes, errno: %d (%s)",lsize + WSIZE,errno,strerror(errno));
#endif
      return DATAIO_ERROR;
    }

  *buffer = buf;
  return data_read(fd,*buffer,MAXBUF);
}

int data_ff(int fd,off_t skipbytes)
{
  file_info_t *infd = getInfo(fd);
  if (infd && infd->sieve) {
    if (infd->writing) {
      fprintf(stderr,"dataIO.data_ff: Error - cannot fast-forward on stream open for writing.\n");
#ifdef USE_ALARM
      esend(ALARM_ERROR,"dataIO.data_ff: cannot fast-forward on stream open for writing.\n");
#endif
      return DATAIO_ERROR;
    }
    else {
      off_t dropbytes = 0;
      while (infd->leftcount <= skipbytes) {
        off_t s;
	int ret;
        char str[80];
        dropbytes += infd->leftcount;
	skipbytes -= infd->leftcount;
	ret = sscanf(fgets(str,80,infd->fp),",%lld\n",&s);
        if (ret == 0 || ret == EOF) {
	  data_flush(fd);
	  sscanf(fgets(str,80,infd->fp),"\n");
          ret = data_read(fd,0,0);
	  if (infd->leftcount == 0) {
	    return DATAIO_EOF;
	  }
	  dropbytes = 0;
	}
	else if (ret == 1) {
          dropbytes += s;
	  if (sscanf(fgets(str,80,infd->fp),"%lld",&infd->leftcount) != 1) {
            fprintf(stderr,"dataIO.data_ff: Error - bad sieve data.\n");
#ifdef USE_ALARM
            esend(ALARM_ERROR,"dataIO.data_ff: bad sieve data.\n");
#endif
	    return DATAIO_ERROR;
	  }
	}
	else {
          fprintf(stderr,"dataIO.data_ff: Error - bad sieve data.\n");
#ifdef USE_ALARM
          esend(ALARM_ERROR,"dataIO.data_ff: bad sieve data.\n");
#endif
	  return DATAIO_ERROR;
	}
      }
      infd->leftcount -= skipbytes;
      return data_ff(infd->link,skipbytes+dropbytes);
    }
  }
  return tape_ff(fd,skipbytes);
}

int data_flush(int fd)
{
  file_info_t *infd = getInfo(fd);
  if (infd && infd->sieve) {
    if (infd->writing) {
      fprintf(infd->fp,"%lld\n",infd->leftcount);
      infd->bytecount += infd->leftcount;
      infd->lastcount = infd->leftcount;
    }
    else if (infd->link) {
      file_info_t *infdin = getInfo(infd->link);
      data_flush(infd->link);
      if (infdin->fp) {
        fclose(infdin->fp);
      }
      file_info[infd->link] = 0;
      free(infdin);
    }
    infd->leftcount = 0;
    infd->link = 0;
    return 0;
  }
  return tape_flush(fd);
}

int data_writeFM(int fd)    /* not to be called on a sieve file! */
{
  tape_flush(fd);
  return tape_writeFM(fd);
}

int data_findFM(int fd,int count)  /* not to be called on a sieve file! */
{
  tape_flush(fd);
  return tape_findFM(fd,count);
}

/* end file */
