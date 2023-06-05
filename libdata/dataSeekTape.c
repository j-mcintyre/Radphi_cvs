/*
  $Log: dataSeekTape.c,v $
  Revision 1.1  1997/05/04 03:21:29  radphi
  Initial revision by lfcrob@dustbunny

*/

static const char rcsid[] = "$Id: dataSeekTape.c,v 1.1 1997/05/04 03:21:29 radphi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>

#include <disData.h>
#include <dataIO.h>
#include <tapeData.h>

static int tapeNo = -1;

static char buf[102400];
static int bufsize = sizeof(buf);

#define FOUND    1
#define NOTFOUND 2
#define MORE     3

static int enableSmartSeeking = 1;

static int doEvent(int ifd,int seekRunNo,const itape_header_t*ih,int *ret_tapeNo)
{
  static int cspill = -1;
  static int crun = -1;
  static int nspills = 0;
  static int lastRun = 0;
  static int repeat = 0;
  int thisRun;

  switch (ih->type)
    {
    default:
      fprintf(stderr,"data.SeekTape: Error: dataIO: Unexpected record type 0x%x ignored\n",ih->type);
      break;

    case TYPE_TAPEHEADER:
      {
	const tapeHeader_t *hdr = (const void*)ih;
	int iTapeNo = (int)hdr->tapeSeqNumber;

	fprintf(stderr,"data.SeekTape: Tape label: Tape %04d written %s",iTapeNo,ctime((const void*)&hdr->startTime));

	if (iTapeNo > 0)
	  {
	    tapeNo = iTapeNo;

	    if (ret_tapeNo) *ret_tapeNo = iTapeNo;
	  }

	break;
      }

    case TYPE_ITAPE:

      fprintf(stderr,"data.SeekTape: Run %4d, spill %d, event %d\n",ih->runNo,ih->spillNo,ih->eventNo);

      thisRun = abs(ih->runNo);

      if (thisRun == lastRun)
	repeat ++;
      else
	repeat = 0;

      if (repeat > 5)
	{
	  fprintf(stderr,"data.SeekTape: Error: run number does not change, give-up\n");
	  return NOTFOUND;
	}

      lastRun = thisRun;
	  
      if (thisRun == seekRunNo)
	{
	  fprintf(stderr,"data.SeekTape: Found: tape %d, run %d\n",tapeNo,seekRunNo);
	  
	  data_findFM(ifd,-1);

	  data_read(ifd,buf,bufsize);

	  return FOUND;
	}

      if ((thisRun < seekRunNo)||(enableSmartSeeking==0))
	{
	  int ret;
	  int runsToSkip = 1;

	  if (enableSmartSeeking)
	    {
	      if (seekRunNo - thisRun > 2)
		{
		  runsToSkip = (seekRunNo - thisRun)*0.66;
		}
	      
	      if (runsToSkip < 1) runsToSkip = 1;
	      if (runsToSkip > 7) runsToSkip = 7;
	    }
	      
	  fprintf(stderr,"data.SeekTape: Seeking %d filemarks forward...\n",runsToSkip);

	  ret = data_findFM(ifd,runsToSkip);

	  if (ret)
	    {
	      if (errno==ENOSPC)
		{
		  if (runsToSkip == 1)
		    {
		      fprintf(stderr,"data.SeekTape: Error: Tape is at EOD (end of data): No more data on this tape\n");
		      return NOTFOUND;
		    }

		  /* we went to far... go back. */

		  fprintf(stderr,"data.SeekTape: We went too far forward, Seeking backward...\n");

		  ret = data_findFM(ifd,-4);

		  enableSmartSeeking = 0;
		}
	    }
	}
      else
	{
	  int ret;

	  fprintf(stderr,"data.SeekTape: Seeking backward...\n");

	  ret = data_findFM(ifd,-4);

	  if (ret)
	    {
	      if (errno==ENOSPC)
		{
		  fprintf(stderr,"data.SeekTape: Error: Tape is at BOT (beginning of tape): No more data on this tape\n");
		  return NOTFOUND;
		}
	    }

	  enableSmartSeeking = 0;
	}

      break;
    }

  return MORE;
}

int data_SeekTape(int ifd,int runNo,int *ret_tapeNo)
{
  itape_header_t *itape = (void*)buf;
  int suspectEOD = 0;
  int ret;

  if (runNo > 0)
    enableSmartSeeking = 1;
  else
    enableSmartSeeking = 0;

  runNo = abs(runNo);

  fprintf(stderr,"data.SeekTape: Looking for run %d\n",runNo);

  while (1)
    {
      ret = data_read(ifd,buf,bufsize);

      if (ret!=DATAIO_EOF) suspectEOD = 0;
      switch (ret)
	{
	default:
	  fprintf(stderr,"data.SeekTape: Error: Unexpected read error code %d\n",ret);
	  return -1;
	case DATAIO_OK:
	  ret = doEvent(ifd,runNo,itape,ret_tapeNo);

	  switch (ret)
	    {
	    case FOUND:    return 0;
	    case NOTFOUND: return -1;
	    case MORE:     break;
	    default:       break;
	    }

	  break;
	case DATAIO_EOF:
	  if (suspectEOD > 5)
	    {
	      fprintf(stderr,"data.SeekTape: Error: Too many Filemarks, assume this is EOD (end of data) and give-up\n");
	      return -1;
	    }
	  suspectEOD += 1;
	  fprintf(stderr,"data.SeekTape: Filemark\n");
	  break;
	case DATAIO_EOT:
	  fprintf(stderr,"data.SeekTape: EOD\n");
	  return -1;
	case DATAIO_BADCRC:
	  fprintf(stderr,"data.SeekTape: Bad CRC\n");
	  return -1;
	case DATAIO_ERROR:
	  fprintf(stderr,"data.SeekTape: Read error\n");
	  return -1;
	}
    }
}

/* end file */
