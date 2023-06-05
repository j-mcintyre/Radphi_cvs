/*
 * evt_dataIO.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <errno.h>

#include <disData.h>
#include <dataIO.h>
#include <tapeIO.h>
#include <crc.h>
#include <ntypes.h>

#include <cefIO.h>
#include <wrapper.h>

#ifdef USE_ALARM
#include <alarm.h>
#endif

#define BUFSIZE 100000 /* Really really big */
#define STREAM_BUFFER_SIZE 524288

#define MY_FOPEN_MAX 128

typedef struct {
   FILE* stream;
   char *streamBuffer;
   int wrongendian;
   intptr_t handle;
} file_info_t;

file_info_t descriptor[MY_FOPEN_MAX];

int evt_open(char *pathname,char *mode)
{
  FILE *fp;
  int fd;
  file_info_t *fdesc;
  unsigned short header[16];
  int headlen=32;

  fp = fopen(pathname,mode);
  if (fp == NULL) {
    perror("Could not open input file");
    fprintf(stderr,"file pathname was: \"%s\"\n",pathname);
    return 0;
  }
  fd = fileno(fp);
  if (fd < 0) {
    fprintf(stderr,"Descriptor to input file \"%s\" is not valid!\n",pathname);
    fclose(fp);
    return 0;
  }
  else if (fd >= MY_FOPEN_MAX) {
    fprintf(stderr,"Maximum open file count %d exceeded!\n",MY_FOPEN_MAX);
    fclose(fp);
    return 0;
  }
  else {
    fdesc = &descriptor[fd];
    fdesc->stream = fp;
    fdesc->streamBuffer = 0;
    fdesc->wrongendian = 0;
    fdesc->handle = 0;
  }	  
  if (*mode == 'r') {
    if (fread(header,1,headlen,fp) != headlen) {
      perror("Could not read from input file");
      fprintf(stderr,"file pathname was: \"%s\"\n",pathname);
      fclose(fp);
      return 0;
    }
    if (header[14] == 0xc0da) {
      if (evOpen(pathname,mode,&fdesc->handle) != S_SUCCESS) {
        perror("Could not open CODA input file");
        fprintf(stderr,"file pathname was: \"%s\"\n",pathname);
        fclose(fp);
        return 0;
      }
    }
    else if (header[14] == 0xdac0) {
      fdesc->wrongendian = 1;
      fprintf(stderr,"evt_dataIO.open: NOTICE: file is wrong-endian, converting on-the-fly.\n");
      if (evOpen(pathname,mode,&fdesc->handle) != S_SUCCESS) {
        perror("Could not open CODA input file");
        fprintf(stderr,"file pathname was: \"%s\"\n",pathname);
        fclose(fp);
        return 0;
      }
    }
    else {
      char version[10];
      char strbuf[30];
      rewind(fp);
      fread(strbuf,4,7,fp);
      if (sscanf(strbuf,"DATA SIEVE version %s",&version) == 1) {
        data_setSieve(fp,0,0);
      }
      rewind(fp);
    }
  }

  if ((fdesc->streamBuffer = malloc(STREAM_BUFFER_SIZE)) == NULL) {
     fprintf(stderr,"evt_open: malloc error creating streamBuffer!\n");
     fclose(fp);
     return 0;
  }
#if 0
  if (setvbuf(fp,fdesc->streamBuffer,_IOFBF,STREAM_BUFFER_SIZE)) {
     fprintf(stderr,"evt_open: setvbuf error opening fd %d!\n",fd);
     fclose(fp);
     return 0;
  }
#endif
  if (fdesc->handle == 0)
     tape_setStream(fp);

  return fd;
}

int evt_data_read(int fd,void *buffer,int bufsize)
{
  int cefBuf[BUFSIZE];
  int retcode;
  file_info_t *fdesc = &descriptor[fd];
  if (fdesc->handle) {
    retcode = evRead(fdesc->handle,cefBuf,BUFSIZE);
    if (retcode == S_EVFILE_UNXPTDEOF) {
      return DATAIO_EOF;
    }
    else if (retcode != S_SUCCESS) {
      return DATAIO_ERROR;
    }

/* When reading a wrongEndian file, a 32-bit byte-swap was already done
 * in evRead.  Most of the words in the raw data record are uint32 so
 * this is OK, with the exception of some short and char members that
 * appear in the CODA bank header structure.  These are patched here.
 * -rtj 29-Mar-2000
 */
    if (fdesc->wrongendian) {
      bankHeader_t *bh;
      endian_convertType1_uint32(&cefBuf[1],sizeof(uint32));
      endian_convertType1_uint16(&cefBuf[1],sizeof(uint16));
      bh = (bankHeader_t *)&cefBuf[1];
      if (bh->dtype == 0x10) {
        endian_convertType1_uint32(&cefBuf[3],sizeof(uint32));
        endian_convertType1_uint16(&cefBuf[3],sizeof(uint16));
      }
    }
    if (wrapEvent(cefBuf,buffer,bufsize)) {
      return DATAIO_ERROR;
    }
    data_clean(buffer);
    return DATAIO_OK;
  }
  else {
    return data_read(fd,buffer,bufsize);
  }
}

int evt_close(int fd)
{
  file_info_t *fdesc = &descriptor[fd];
  if (fdesc->handle) {
    evClose(fdesc->handle);
  }
  else {
    data_flush(fd);
  }
  fclose(fdesc->stream);
  free(fdesc->streamBuffer);
  return 0;
}

/* end file */
