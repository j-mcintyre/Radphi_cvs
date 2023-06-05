/*
 * disIO.c
*/

static const char sccsid[] = "@(#)"__FILE__"\t\t4.4\tCreated 4/22/97 17:41:42, \tcompiled "__DATE__;

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <netdb.h>

#include <disIO.h>
#include <communications.h>

#include <signal.h>

#ifdef IRIX

#include <bstring.h>

#else

#if (!defined(HPUX))
#include <sys/select.h>
#endif

#endif

#ifdef HPUX
extern int h_errno;
#endif

#ifdef SOLARIS
extern int h_errno;
#endif

int disIO_socket = (-1);

uint32 disIO_getaddress(const char*server)
{
    struct hostent *host;
    uint32 netaddress;    /* address in network format */
    uint32 address = 0;   /* address in local format */
    
    if (server) host = gethostbyname(server);
    else host = gethostent();
    if (host) 
      {
	netaddress = *(uint32*)*(host->h_addr_list);
	address = ntohl(netaddress);

#ifdef UNDEF
	fprintf(stderr,"disIO.getaddress: Lookup for host (%s) gives %d.%d.%d.%d\n",server,
		(address>>24)&0xff,(address>>16)&0xff,
		(address>>8)&0xff,(address)&0xff);
#endif
      }
    else {
	fprintf(stderr,"disIO_getaddress: Error: Cannot resolve hostname '%s': h_errno=%d\n",server,h_errno);
	return 0;
      }
    return(address);
}

const char* disIO_getHost(const char*server)
{
  int i;
  static char *host = NULL;
  uint32 hostLength = 0;

  if (strlen(server) >= hostLength)
    {
      if (host) free(host);
      hostLength = strlen(server) + 1;
      host = malloc(hostLength);
    }

  for (i=0; server[i]!=0; i++)
    {
      if (server[i] == ':') break;
      host[i] = server[i];
    }

  host[i] = 0;

  return host;
}

int disIO_getPort(const char*server)
{
  int i;

  for (i=0; server[i]!=0; i++) if (server[i] == ':') break;

  if (server[i] != ':') return -1;

  i++;

  return (int)strtoul(&server[i],NULL,0);
}

#define TCP 6

#define MINTIME  5 /* seconds */
#define MAXTIME 30 /* seconds */

/*static void (*saved_handler)(int v);*/

int disIO_connect(const char*server,int defaultHandle)
{
  const char* host;
  int handle = 0;
  int ms;
  unsigned int IPaddr = 0;

  if (disIO_socket >= 0) disIO_disconnect();

  if (server == NULL)
    {
      fprintf(stderr,"disIO_connect: The server name is NULL\n");
      return -1;
    }

  signal(SIGPIPE,SIG_IGN);

  host   = disIO_getHost(server);
  handle = disIO_getPort(server);
  IPaddr = disIO_getaddress(host);

  if (handle < 0) handle = defaultHandle;
  
  fprintf(stderr,"disIO.connect: Connecting to the Dispatcher on %s (%d.%d.%d.%d), TCP port %d...",
	  host,
	  (IPaddr>>24)&0xff,
	  (IPaddr>>16)&0xff,
	  (IPaddr>>8)&0xff,
	  (IPaddr)&0xff,
	  handle);
  
  fflush(stderr);

  {
    struct protoent  *protocol;

    protocol = getprotobyname("tcp");
    if (protocol == NULL)
      {
	fprintf(stderr,"disIO.connect: getprotobyname('tcp') failed, errno: %d (%s)\n",errno,strerror(errno));
	return -1;
      }

    ms = socket(PF_INET,SOCK_STREAM,protocol->p_proto);
    if (ms<0)
      {
	fprintf(stderr,"disIO.connect: socket(PF_INET,SOCK_STREAM,tcp) failed, errno: %d (%s)\n",errno,strerror(errno));
	return -1;
      }
  }

  {
    int ret;
    struct sockaddr_in in_name;

    if (IPaddr == 0)
      {
	close(ms);
	return -1;
      }

    in_name.sin_family = AF_INET;
    in_name.sin_port = htons(handle);
    in_name.sin_addr.s_addr = (int)htonl(IPaddr);
  
    ret = connect(ms,(void*)&in_name,sizeof(in_name));
    if (ret<0)
      {
	fprintf(stderr,"disIO.connect: connect to %s:%d failed, errno: %d (%s)\n",host,handle,errno,strerror(errno));
	close(ms);
	return -1;
      }
  }

  {
    int err;
#ifdef HPUX
    int soBufSize = 58254;
#else
    int soBufSize = 65000;
#endif

    err = setsockopt(ms,SOL_SOCKET,SO_SNDBUF,(const char*)&soBufSize,sizeof(soBufSize));
    if (err == (-1))
      {
	fprintf(stderr,"disIO.connect: setsockopt(SO_SNDBUF) failed, errno: %d (%s)\n",errno,strerror(errno));
      }
    
    err = setsockopt(ms,SOL_SOCKET,SO_RCVBUF,(const char*)&soBufSize,sizeof(soBufSize));
    if (err == (-1))
      {
	fprintf(stderr,"disIO.connect: setsockopt(SO_RCVBUF) failed, errno: %d (%s)",errno,strerror(errno));
      }
  }
    
  fprintf(stderr,"connected.\n");

  disIO_socket = ms;
  return 0;
}

int disIO_disconnect(void)
{
  if (disIO_socket>=0)
    {
      close(disIO_socket);
      disIO_socket = (-1);
    }
  return 0;
}

static int protected_send(int fd,const void*ptr,int size)
{
  int ret;
  int32 wr;

  ret = (-1);
  wr = write(fd,ptr,size);
  if (wr == size) ret = 0;
  else if (wr<0)
    {
      if (errno == EPIPE)        ret = DISIO_EOF;
      if (errno == ETIMEDOUT)    ret = DISIO_EOF;
      if (errno == ECONNRESET)   ret = DISIO_EOF;
      if (errno == ECONNREFUSED) ret = DISIO_EOF;
    }
  else if (wr != size) 
    { 
      fprintf(stderr,"disIO_wr: Short send(fd=%d,ptr=0x%x,size=%d), errno: %d (%s)\n",fd,(intptr_t)ptr,size,errno,strerror(errno)); 
      ret = (-1); 
    }

  return ret;
}

int disIO_sendPacket(int Wsocket,int type_,const void *ptr,int size)
{
  int ret;
  DatagramHeader_t header;

  header.nbytes = htonl(size);
  header.transcomputer = htonl(TransComputerCode);
  header.type = htonl(type_);

  ret = protected_send(disIO_socket,&header,sizeof(header));
  if (ret!=0) { 
    fprintf(stderr,"disIO_sendPacket: Header send() failed: Wsocket: %d, type: %d, ptr: 0x%x, size: %d, errno: %d (%s)\n",
	    Wsocket,type_,(intptr_t)ptr,size,errno,strerror(errno)); 
    return ret; 
  }

  ret = protected_send(disIO_socket,ptr,size);
  if (ret<0) { 
    fprintf(stderr,"disIO_sendPacket: Data send() failed: Wsocket: %d, type: %d, ptr: 0x%x, size: %d, errno: %d (%s)\n",
	    Wsocket,type_,(intptr_t)ptr,size,errno,strerror(errno)); 
    return ret;
  }

  return 0;
}

int disIO_cmd(int Wsocket, const char*CommandText)
{
  return disIO_sendPacket(Wsocket,COMMAND,CommandText,(int)strlen(CommandText)+1);
}

int disIO_command(const char*CommandText)
{
  return disIO_cmd(disIO_socket,CommandText);
}

int disIO_checkSocket(int fd)
{
  struct timeval timeint;
  fd_set readset;
  int nfound;

  FD_ZERO(&readset);
  FD_SET(fd,&readset);

  timeint.tv_sec = 0;
  timeint.tv_usec = 0;

  nfound = select(fd + 1,&readset,0,0,&timeint);
  if (nfound<0) { 
    perror("disIO_checkSocket: select() failed"); 
    return -1; 
  }

  return FD_ISSET(fd,&readset);
}

static int read1(int fd,void*ptr_,int count)
{
  char*ptr = (char*)ptr_;
  int32 nread;
  int n=0;

  while (count)
    {
      nread = read(fd,ptr,count);
      if (nread<0) 
	{ 
	  fprintf(stderr,"disIO_rd: read(fd=%d,ptr=0x%x,count=%d) failed, errno: %d (%s)\n",
		  fd,(intptr_t)ptr,count,errno,strerror(errno)); 
	  if (errno==EPIPE)        return 0;
	  if (errno==ECONNREFUSED) return 0;
	  if (errno==ECONNRESET)   return 0;
	  if (errno==ETIMEDOUT)    return 0;
	  if (errno==EINTR)        continue;
	  return -1; 
	}
      if (nread==0) return n; 

      ptr += nread;
      count -= nread;
      n += nread;
    }

  return n;
}

int disIO_readHeader(int fd,DatagramHeader_t *header,int wait)
{
  int ret;
  if (wait==0)
    {
      fd_set readset;
      struct timeval timeint;
      
      FD_ZERO(&readset);
      FD_SET(fd,&readset);

      timeint.tv_sec = 0;
      timeint.tv_usec = 0;

      ret = select(fd + 1,&readset,0,0,&timeint);
      if (ret==0) return 1;
      if (ret<0) { 
	perror("disIO_rd: select() failed");
	return -1; 
      }
    }

  ret = read1(fd,header,sizeof(*header));
  if (ret==0) return DISIO_EOF;
  if (ret!=sizeof(*header)) { 
    fprintf(stderr,"disIO_rd: Bad header read(): See prev message\n"); 
    return -1; 
  }

/*
  fprintf(stderr,"disIO_readHeader: %d bytes %d type\n",ntohl(header->nbytes),ntohl(header->type));
*/
  return 0;
}

int disIO_readPacket(int Rsocket,void**buffer,int*bufsize,int wait)
{
  DatagramHeader_t header;
  int size,type;
  int ret;

  ret = disIO_readHeader(Rsocket,&header,wait);
  if (ret<0)
    {
      fprintf(stderr,"disIO_readPacket: Header read failed: See prev message\n"); 
      return -1;
    }
  if (ret==1) return DISIO_NODATA;
  if (ret==DISIO_EOF) return ret;

  size = ntohl(header.nbytes);
  type = ntohl(header.type);

  if (*buffer==NULL)
    {
      *buffer = malloc(size);
      if (*buffer==NULL)
	{
	  perror("disIO_readPacket: malloc() failed");
	  fprintf(stderr,"disIO_readPacket: Cannot allocate %d bytes for data buffer: See prev message\n",size);
	  return -1;
	}
      *bufsize = size;
    }
  else
    {
      if (size>*bufsize)
	{ 
	  fprintf(stderr,"disIO_readPacket: Data (%d bytes) longer than buffer (%d bytes)\n",size,*bufsize);
	  return DISIO_OVERRUN;
	}
    }

  ret = read1(Rsocket,*buffer,size);
  if (ret<0) 
    { 
      fprintf(stderr,"disIO_readPacket: Data read(%d bytes) failed: See previous messages",size); 
      return -1; 
    }
  if (ret==0)
    {
      fprintf(stderr,"disIO_readPacket: Data read(%d bytes): zero read, See previous messages",size); 
      return DISIO_EOF;
    }

  return type;
}

int disIO_readRAW(void*buffer,int bufsize,int *nreadPtr,int wait)
{
  int Rsocket = disIO_socket;
  int ret;
  int nread = 0;

  if (nreadPtr == NULL)
    nreadPtr = &nread;

  *nreadPtr = bufsize;
  ret = disIO_readPacket(Rsocket,&buffer,nreadPtr,wait);
  if (ret<0) return ret;
  if (ret == DISIO_NODATA) return ret;
  if (ret == DISIO_EOF) return ret;
  if (ret == COMMAND) return DISIO_COMMAND;
  if (ret != DATA_RECORD) return DISIO_WRONGTYPE;

  return DISIO_OK;
}

int disIO_readRAW_alloc(void**buffer,int *nreadPtr,int wait)
{
  int Rsocket = disIO_socket;
  int ret;
  int nread = 0;

  if (nreadPtr == NULL)
    nreadPtr = &nread;

  *buffer = NULL;
  *nreadPtr = 0;
  ret = disIO_readPacket(Rsocket,buffer,nreadPtr,wait);
  if (ret<0) return ret;
  if (ret == DISIO_NODATA) return ret;
  if (ret == DISIO_EOF) return ret;
  if (ret == COMMAND) return DISIO_COMMAND;
  if (ret != DATA_RECORD) return DISIO_WRONGTYPE;

  return DISIO_OK;
}

int disIO_readCOOKED(void*buffer,int bufsize,int *nreadPtr,int wait)
{
  int Rsocket = disIO_socket;
  int ret;
  int nread = 0;

  if (nreadPtr == NULL)
    nreadPtr = &nread;

  *nreadPtr = bufsize;
  ret = disIO_readPacket(Rsocket,&buffer,nreadPtr,wait);
  if (ret<0) return ret;
  if (ret == DISIO_NODATA) return ret;
  if (ret == DISIO_EOF) return ret;
  if (ret == COMMAND) return DISIO_COMMAND;
  if (ret != ITAPE_RECORD) return DISIO_WRONGTYPE;
  return DISIO_OK;
}

int disIO_readCOOKED_alloc(void**buffer,int *nreadPtr,int wait)
{
  int Rsocket = disIO_socket;
  int ret;
  int nread = 0;

  if (nreadPtr == NULL)
    nreadPtr = &nread;

  *buffer = NULL;
  *nreadPtr = 0;
  ret = disIO_readPacket(Rsocket,buffer,nreadPtr,wait);
  if (ret<0) return ret;
  if (ret == DISIO_NODATA) return ret;
  if (ret == DISIO_EOF) return ret;
  if (ret == COMMAND) return DISIO_COMMAND;
  if (ret != ITAPE_RECORD) return DISIO_WRONGTYPE;
  return DISIO_OK;
}

int disIO_writeRAW(const void*data,int size)
{
  return disIO_sendPacket(disIO_socket,DATA_RECORD,data,size);
}

int disIO_writeCOOKED(const void*data,int size)
{
  return disIO_sendPacket(disIO_socket,ITAPE,data,size);
}

int disIO_read_fd(void)
{
  return disIO_socket;
}


int disIO_mapRequest(int mapRequest,   /* from the MAP_REQUEST_XXX set */
		     const char *systemName,
		     const char *subsystemName,
		     const char *itemName,
		     int accessKey,
		     int dataType,     /* from the MAP_DATATYPE_XXX set */
		     int dataCount,
		     const void*dataPointer)
{
  int ret;
  uint32 reqSize = sizeof(map_request_t)
    + strlen(systemName) + 1
    + strlen(subsystemName) + 1
    + strlen(itemName) + 1
    + sizeof(float) * dataCount;

  map_request_t *req = malloc(reqSize);

  uint32 systemName_off,data_off;
  uint32 subsystemName_off,itemName_off;
  uint32 dataSz,dataEnd;
  int reqSz;

  req->requestCode = htonl(mapRequest);
  req->dataType    = htonl(dataType);
  req->dataCount   = htonl(dataCount);
  req->accessKey   = htonl(accessKey);

  systemName_off    = 0;
  subsystemName_off = systemName_off    + strlen(systemName)    + 1;
  itemName_off      = subsystemName_off + strlen(subsystemName) + 1;
  data_off          = itemName_off      + strlen(itemName)      + 1;

  strcpy(&req->data[systemName_off],    systemName);
  strcpy(&req->data[subsystemName_off], subsystemName);
  strcpy(&req->data[itemName_off],      itemName);

  data_off = ((data_off/4)+1)*4; /* align data on a 4-bytes boundary */

  switch (dataType)
    {
    case MAP_DATATYPE_INT:   dataSz = sizeof(int);   break;
    case MAP_DATATYPE_FLOAT: dataSz = sizeof(float); break;
    case MAP_DATATYPE_STRING: dataSz = sizeof(char); break;
    default:
      fprintf(stderr,"disIO.mapRequest: Error: Unknown data type: %d. Request ignored\n",dataType);
      return DISIO_WRONGTYPE;
    }

  memcpy(&req->data[data_off],dataPointer,dataCount * dataSz);

  dataEnd = data_off + dataCount * dataSz;

  req->systemName_offset    = (int)htonl(systemName_off);
  req->subsystemName_offset = (int)htonl(subsystemName_off);
  req->itemName_offset      = (int)htonl(itemName_off);
  req->data_offset          = (int)htonl(data_off);

  reqSz = (int) ((char*)&req->data[dataEnd] - (char*)req);

  ret = disIO_sendPacket(disIO_socket,MAP_REQUEST,req,reqSz);
  free(req);
  return ret;
}

/* end file */
