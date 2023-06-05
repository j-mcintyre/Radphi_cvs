/*
  $Log: disIO.h,v $
  Revision 1.1  1997/05/04 03:20:25  radphi
  Initial revision by lfcrob@dustbunny
  Dispatcher I/O library

*/

#include <ntypes.h>

#ifndef _DISIO1_INCLUDED
#define _DISIO1_INCLUDED

#define DISIO_OK 0
#define DISIO_EOF     100   /* socket connection was closed from the other side */
#define DISIO_NODATA  101   /* no data is available from the dispatcher yet */
#define DISIO_COMMAND 102   /* received a command */
#define DISIO_OVERRUN 103   /* supplied buffer too small for data (fatal error) */
#define DISIO_WRONGTYPE 104 /* the Dispatcher sent us data of an unexpected type */

#define DISIO_Handle_RADPHI 10357

extern int disIO_socket; /* the default Dispatcher connection socket */

/*
 * if 'server' is in the 'host:port' format, the 'port' argument is not used 
*/ 

int disIO_connect(const char*server,int port);
int disIO_disconnect(void);

/*
 * see the disIO manual for a list of commands
*/

int disIO_command(const char*str);

/*
 * Low level communication routines.
 * The 'socket' argument specifies the Dispatcher connection to use.
*/

uint32 disIO_getaddress(const char*hostName); /* returns the IP address in host format or 0 if error */

int disIO_sendPacket(int Wsocket,int packetType,const void *ptr,int size);
int disIO_readPacket(int Rsocket,void**buffer,int*bufsize,int wait);
int disIO_cmd(int Wsocket, const char*CommandText);
int disIO_checkSocket(int Rsocket); /* check if data is available. Return value: 1 if yes, 0 if no and -1 if error */

const char* disIO_getHost(const char*server); /* get the 'host' part of the server name */
        int disIO_getPort(const char*server); /* get the 'port' part of the server name */
/*
 * High level routines
*/

int disIO_read_fd(void); /* get the Dispatcher connection socket */
int disIO_readRAW(void*buffer,int bufsize,int *nread,int wait);
int disIO_readCOOKED(void*buffer,int bufsize,int *nread,int wait);

int disIO_readRAW_alloc(void**buffer,int *bufsize,int wait);
int disIO_readCOOKED_alloc(void**buffer,int *bufsize,int wait);

int disIO_writeRAW(const void*buffer,int bufsize);
int disIO_writeCOOKED(const void*buffer,int bufsize);

/*
 * This implements a set of Map Manager 'write-only' requests.
 *
 * The 'systemName', 'subsystemName', 'itemName' and 'accessKey'
 * correspond to the 'file', 'subsystem', 'item' and 'firsttime'
 * arguments in the map_put_XXX(...) calls.
*/

#define MAP_REQUEST_PUT  102
#define MAP_REQUEST_ADD  103
#define MAP_REQUEST_ZERO 104

#define MAP_DATATYPE_INT    201
#define MAP_DATATYPE_FLOAT  202
#define MAP_DATATYPE_STRING 203

int disIO_mapRequest(int mapRequest,   /* from the MAP_REQUEST_XXX set */
		     const char *systemName,
		     const char *subsystemName,
		     const char *itemName,
		     int accessKey,
		     int dataType,     /* from the MAP_DATATYPE_XXX set */
		     int dataCount,
		     const void*dataPointer);

#endif
/* end file */
