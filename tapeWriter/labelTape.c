/*
*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#include <disData.h>
#include <dataIO.h>
#include <tapeIO.h>
#include <tapeData.h>

#include <mtime.h>

void usage(void)
{
  fprintf(stderr,"Usage: tapeNo [-tXXXX] [-e] [-r] tapedevice\n");
  fprintf(stderr,"Switches:\n");
  fprintf(stderr,"\t-tXXXX  is the tape number\n");
  fprintf(stderr,"\t-e      eject the tape after the operation\n");
  fprintf(stderr,"\t-r      rewind the tape before the operation\n");
  fprintf(stderr,"\t-x      test the tape quality before writing the tape label\n");
  exit(1);
}

int debug=0;
int exitFlag = 0;

int main(int argc,char*argv[])
{
  int ifd;
  int i;
  char*arg;
  int tapeNo = 0;
  int autoEject = 0;
  int autoRewind = 0;
  int testTape = 0;

  if (argc < 2) usage();

  for (i=1; i<argc; i++)
    {
      if ((*argv[i])=='-')
	{
	  arg = &argv[i][2];
	  switch (argv[i][1])
	    {
	    case 'e': autoEject = 1; break;
	    case 'r': autoRewind = 1; break;
	    case 'x': testTape = 1; break;
	    case 't': tapeNo = strtol(arg,NULL,0); break;
	    case 'h': usage(); break;
	    default:
	      fprintf(stderr,"Unknown option %s\n",argv[i]);
	      usage();
	      break;
	    }
	}
    }

  for (i=1; i<argc; i++)
    if (*argv[i]!='-')
      {

	if (tapeNo == 0)
	  {
	    ifd = open(argv[i],O_RDONLY);
	  }
	else
	  {
	    ifd = open(argv[i],O_RDWR);
	  }

	if (ifd < 0)
	  {
	    fprintf(stderr,"tapeNo: Cannot open tape device %s: %s\n",argv[i],strerror(errno));
	    exit(1);
	  }

	break;
      }

  if (tapeNo == 0)  /* just read the current tape label */
    {
      itape_header_t *itape = NULL;
      int ret;

      if (autoRewind) tape_rewind(ifd);

      ret = data_read_alloc(ifd,(void*)&itape);

      switch (ret)
	{
	default:
	  fprintf(stderr,"readEvent: Unknown data_read_alloc() return code %d\n",ret);
	  exitFlag = 1;
	  break;
	case DATAIO_OK:

	  switch (itape->type)
	    {
	    default:
	      fprintf(stderr,"dataIO: Unexpected rec type 0x%x ignored\n",itape->type);
	      break;

	    case TYPE_TAPEHEADER:
	      {
		const tapeHeader_t *hdr = (const void*)itape;
		int iTapeNo = hdr->tapeSeqNumber;

		printf("Tape %04d, label written %s",iTapeNo,ctime((const void*)&hdr->startTime));

		break;
	      }

	    case TYPE_ITAPE:
	      fprintf(stderr,"readEvent: No tape header\n");
	      break;

	    }
	  break;
	case DATAIO_EOF:
	  fprintf(stderr,"readEvent: EOF\n");
	  break;
	case DATAIO_BADCRC:
	  fprintf(stderr,"readEvent: Bad CRC\n");
	  break;
	case DATAIO_ERROR:
	  fprintf(stderr,"readEvent: dataIO.read_alloc() failed\n");
	  break;
	}
    }
  else /* write a new tape label */
    {
      int tapeIsOk = 1;

      if (testTape)
	{
	  char buf[102400];
	  int bufsize = sizeof(buf);
	  int wr = 0;
	  unsigned long start,endt;
	  int i;
	  int count = 200;
	  double kps = 0;

	  errno = 0;

	  fprintf(stderr,"tapeNo: Testing the tape, please wait.\n");
	  fprintf(stderr,"tapeNo: Loading...");

	  start = (int)mtime();
	  wr = write(ifd,buf,bufsize);

	  if (wr != bufsize)
	    {
	      fprintf(stderr,"tapeNo: Error: write(%d) returned %d, errno: %s\n",
		      bufsize,wr,strerror(errno));
	      exit(1);
	    }

	  endt = mtime();

	  fprintf(stderr,"%d msec\n",endt - start);
	  fprintf(stderr,"tapeNo: Writing %d data blocks...",count);

	  start = mtime();

	  for (i=1; i<=count; i++)
	    {
	      wr = write(ifd,buf,bufsize);
	      
	      if (wr != bufsize)
		{
		  fprintf(stderr,"tapeNo: Error: write(%d) returned %d, errno: %s\n",
			  bufsize,wr,strerror(errno));
		  exit(1);
		}

	      if ((i%20)==0)
		fprintf(stderr,".%d",i);
	    }

	  endt = mtime();

	  kps = (bufsize*count)/(double)(endt-start);

	  fprintf(stderr,". %d msec, %6.1lf Kbytes/sec\n",endt - start,(double)kps);
	}

      if (tapeIsOk)
	{
	  tapeHeader_t hdr;
	  
	  hdr.length = sizeof(hdr) - 4;
	  hdr.type = TYPE_TAPEHEADER;
	  hdr.transComputerCode = 1;
	  hdr.ngroups = 1;
	  hdr.tapeSeqNumber = tapeNo;
	  hdr.startTime = time(NULL);
	  sprintf(hdr.devname,"tapeNo");
	  sprintf(hdr.from,"tapeNo");
	  
	  fprintf(stderr,"tapeNo: Rewinding...\n");
	  tape_rewind(ifd);
	  fprintf(stderr,"tapeNo: Writing the label...\n");
	  data_write(ifd,&hdr);
	  data_flush(ifd);

	  if (autoEject)
	    {
	      tape_writeFM(ifd);
	      tape_unload(ifd);
	    }

	  close(ifd);

	  printf("Tape %04d, label written %s",hdr.tapeSeqNumber,ctime((const void*)&hdr.startTime));

	  exit(0);
	}
    }

  if (autoEject) tape_unload(ifd);
  close(ifd);

  exit(0);
}

/* end file */
