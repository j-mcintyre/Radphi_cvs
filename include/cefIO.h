

#include <evfile_msg.h>
#ifndef EFIO_H_INCLUDED
#define EFIO_H_INCLUDED

#include <inttypes.h>

/* $Log: cefIO.h,v $
/* Revision 1.2  1997/01/31 17:28:12  lfcrob
/* Added include of evfile_msg.h
/*
 * Revision 1.1  1997/01/31  16:34:23  lfcrob
 * Initial revision
 *

 */

int evOpen(char *filename,char *flags,intptr_t *handle);
int evRead(intptr_t handle,int *buffer,int buflen);
int evWrite(intptr_t handle,int *buffer);
int evIoctl(intptr_t handle,char *request,void *argp);
int evClose(intptr_t handle);
#endif
