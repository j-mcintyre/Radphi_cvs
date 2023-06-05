/* 
   $Log: wrapFile.c,v $
   Revision 1.5  2014/01/06 04:58:42  radphi
   * Makefile - changes to build better on Centos 6 [rtj]
   * wrapFile.c - changes to get rid of compiler warnings [rtj]

   Revision 1.4  1997/07/01 13:50:03  radphi
   Changes made by lfcrob@jlabs2
   Added stdout

 * Revision 1.3  1997/05/08  04:55:18  radphi
 * Changes made by lfcrob@dustbunny
 * Now reports on different "data" triggers
 *
 * Revision 1.2  1997/03/26  15:20:42  lfcrob
 * Fixed "-r" bug
 *
 * Revision 1.1  1997/02/12  15:15:06  lfcrob
 * Initial revision
 *
   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <eventType.h>
#include <triggerType.h>
#include <disData.h>
#include <wrapper.h>
#include <cefIO.h>
#include <dataIO.h>

static char rcsid[] = "$Id: wrapFile.c,v 1.5 2014/01/06 04:58:42 radphi Exp $";

#define BUFSIZE 100000 /* Really really big */

int main(int argc, char *argv[]);
void Usage(void);

void Usage(void)
{
  fprintf(stderr,"wrapFile <options> cefFile itapeFile\n");
  fprintf(stderr,"Options:\n");
  fprintf(stderr,"\t-r #\t\tSet default run number to #.\n");
  fprintf(stderr,"\t-h\t\tPrint this message.\n");
}


int main(int argc, char *argv[])
{
  int iarg;
  char *argptr;

  char *inputFile=NULL;
  char *outputFile=NULL;
  itape_header_t *itape;
  int *cefBuf;
  intptr_t cefHandle;
  FILE *fpOut;
  int defRun;
  int nProcessed=0;

  int nData=0;
  int nSync=0;
  int nPreStart=0;
  int nGo=0;
  int nPause=0;
  int nEnd=0;
  int nUnknown=0;
  int nTrigData=0;
  int nTrigLGDMon=0;
  int nTrigRPDMon=0;
  int nTrigPed=0;
  int nTrigBT=0;
  int nTrigUnknown=0;

  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'r':
	defRun = atoi(argv[++iarg]);
	wrapperSetDefaults( WRAP_DEF_RUN,defRun);
	break;
      case 'h':
	Usage();
	return(0);
      default:
	fprintf(stderr,"Unknown argument : %s\n",argv[iarg]);
	Usage();
	return(1);
      }
    }
    else{
      if(!inputFile)
	inputFile = argv[iarg];
      else{
	if(!outputFile)
	  outputFile = argv[iarg];
	else{
	  Usage();
	  return(1);
	}	
      }
    }
  if(!inputFile){
    Usage();
    return(1);
  }
  
  /* Ready to rock 'n roll */

  if(evOpen(inputFile,"r",&cefHandle) == S_SUCCESS){
    cefBuf = malloc(BUFSIZE);
    if(outputFile)
      fpOut = fopen(outputFile,"w");
    else
      fpOut = stdout;
    if(outputFile && !fpOut){
      fprintf(stderr,"Unable to open output file %s\n",outputFile);
      exit(1);
    }
    else{
      itape= malloc(BUFSIZE);
      while((evRead(cefHandle,cefBuf,BUFSIZE)  == S_SUCCESS)){
	if(!wrapEvent(cefBuf,itape,BUFSIZE)){
	  data_clean(itape);
	  data_write(fileno(fpOut),itape);
	  switch(itape->eventType){
	  case EV_DATA:
	    switch(itape->trigger){
	    case TRIG_DATA:
	      nTrigData++;
	      break;
	    case TRIG_LGDMON:
	      nTrigLGDMon++;
	      break;
	    case TRIG_PED:
	      nTrigPed++;
	      break;
	    case TRIG_BT:
	      nTrigBT++;
	      break;
	    case TRIG_RPDMON:
	      nTrigRPDMon++;
	      break;
	    default:
	      nTrigUnknown++;
	      break;
	    }
	    nData++;
	    break;
	  case EV_SYNC:
	    nSync++;
	    break;
	  case EV_PRESTART:
	    nPreStart++;
	    break;
	  case EV_GO:
	    nGo++;
	    break;
	  case EV_PAUSE:
	    nPause++;
	    break;
	  case EV_END:
	    nEnd++;
	    break;
	  case EV_UNKNOWN:
	    nUnknown++;
	    break;
	  default:
	    fprintf(stderr,"EVent type : %d\n",itape->eventType);
	  }	  
	  if((++nProcessed%10) == 0){
	    fprintf(stderr,"%d\r",nProcessed);
	    fflush(stderr);
	  }
	}
	else{
	  fprintf(stderr,"Failed to wrap event!\n");
	}
      }
      data_flush(fileno(fpOut));
      fclose(fpOut);
    }
  }
  else{
    fprintf(stderr,"Could not open %s for reading.\n",inputFile);
    return(1);
  }
  fprintf(stderr,"Processed %d events\n",nProcessed);
  fprintf(stderr,"\t%d data\n",nData);
  fprintf(stderr,"\t\t%d physics events\n",nTrigData);
  fprintf(stderr,"\t\t%d LGD monitor\n",nTrigLGDMon);
  fprintf(stderr,"\t\t%d RPD monitor\n",nTrigRPDMon);
  fprintf(stderr,"\t\t%d pedestal\n",nTrigPed);
  fprintf(stderr,"\t\t%d basetest\n",nTrigBT);
  fprintf(stderr,"\t\t%d unknown\n",nTrigUnknown);
  fprintf(stderr,"\t%d sync\n",nSync);
  fprintf(stderr,"\t%d prestart\n",nPreStart);
  fprintf(stderr,"\t%d go\n",nGo);
  fprintf(stderr,"\t%d pause\n",nPause);
  fprintf(stderr,"\t%d end\n",nEnd);
  fprintf(stderr,"\t%d unknown\n",nUnknown);
  return(0);
}










