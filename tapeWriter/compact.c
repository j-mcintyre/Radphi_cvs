/*
 * compact
 *
 * A standalone program to input a series of raw data files and generate
 * on the standard output a compacted stream of itape records.  Special
 * records are copied verbatim, but TRIG_DATA records are unpacked and
 * the GROUP_RAW data group is replaced with a GROUP_FROZEN digest of its
 * contents.  The data volume in the output stream is typically 75% that
 * of the input stream.  Besides reducing the volume, compact translates
 * the data to the native endianness of the local machine, hence improving
 * efficiency if several analysis analysis passes are subsequently made
 * over the * data.  Subsequent analysis that uses the unpackEvent()
 * function to access raw data will work transparently on the compact
 * output stream without modification.
 *
 * Richard Jones, July 10 2000
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <eventType.h>
#include <triggerType.h>
#include <disData.h>
#include <dataIO.h>
#include <makeHits.h>
#include <tapeData.h>
#include <math.h>
#include <unpackData.h>


int main(int argc, char *argv[]);
void Usage(void);

itape_header_t event[BUFSIZE];

void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  compact <options> file1 file2 ...\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-h\t\tPrint this message\n");
  fprintf(stderr,"\t-M#\t\tProcess at most # events.\n");
}

int main(int argc, char *argv[])

{
  int iarg;
  char *argptr;
  int fpIn;
  int maxEvents=0;
  int srun[]={0,0};
  int nEvents=0;
  int nData=0;
  int nTrigData=0;
  int nTrigLGDMon=0;
  int nTrigBGVMon=0;
  int nTrigBT=0;
  int nTrigPed=0;
  int nTrigUnknown=0;
  int nSync=0;
  int nPreStart=0;
  int nGo=0;
  int nPause=0;
  int nEnd=0;
  int nUnknown=0;

  tape_setStream(stdout);

  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'M':
	maxEvents = atoi(++argptr);
	break;
      case 'h':
	Usage();
	exit(1);
      default:
	fprintf(stderr,"Unknown argument: %s\n",argv[iarg]);
	Usage();
	exit(1);
      }
    }
  
  for(iarg=1;iarg<argc;iarg++){
    if(*argv[iarg] != '-'){
      if( (fpIn = evt_open(argv[iarg],"r")) == 0){
	fprintf(stderr,"Failed to open %s for reading\n",argv[iarg]);
	exit(1);
      }
      else{
	fprintf(stderr,"Reading from %s\n",argv[iarg]);
	while( (maxEvents ? maxEvents > nEvents : 1)  &&
               (evt_data_read(fpIn,event,BUFSIZE) == DATAIO_OK) ){
	  if( (++nEvents%1000) == 0){
	    fprintf(stderr,"%d\r",nEvents);
	    fflush(stderr);
	  }
          if (event->runNo != srun[0]) {
            if (setup_makeHits(event)) {
              fprintf(stderr,"setup_makeHits failed, quitting!\n");
              exit(1);
            }
            srun[0] = event->runNo;
            srun[1] = 0;
          }
          switch(event->eventType) {
          case EV_DATA:
            switch(event->trigger) {
            case TRIG_DATA:
              nTrigData++;
              unpackEvent(event,BUFSIZE);
              freezeData(event,BUFSIZE);
              if (data_getGroup(event,GROUP_SCALERS,0)) {
                data_removeGroup(event,GROUP_FROZEN,0);
                data_clean(event);
                event->trigger = 0;
                data_write(1,event);
                data_removeGroup(event,GROUP_SCALERS,0);
                freezeData(event,BUFSIZE);
                event->trigger = TRIG_DATA;
              }
              data_removeGroup(event,GROUP_RAW,0);
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
            case TRIG_BGVMON:
              nTrigBGVMon++;
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
            fprintf(stderr,"EVent type : %d\n",event->eventType);
          }	  
          data_clean(event);
          data_write(1,event);
          srun[1]++;
	}
	if(fpIn){
	  evt_close(fpIn);
	  fpIn = 0;
	}
      }      
    }  
  }  
  fprintf(stderr,"Processed %d records:\n",nEvents);
  fprintf(stderr,"\t%d data",nData);
  fprintf(stderr,"\t%d physics",nTrigData);
  fprintf(stderr,"\t%d LGD monitor",nTrigLGDMon);
  fprintf(stderr,"\t%d BGV monitor\n",nTrigBGVMon);
  fprintf(stderr,"\t\t\t%d pedestal",nTrigPed);
  fprintf(stderr,"\t%d basetest",nTrigBT);
  fprintf(stderr,"\t%d unknown\n",nTrigUnknown);
  fprintf(stderr,"\t%d sync",nSync);
  fprintf(stderr,"\t%d prestart",nPreStart);
  fprintf(stderr,"\t%d go",nGo);
  fprintf(stderr,"\t%d pause",nPause);
  fprintf(stderr,"\t%d end",nEnd);
  fprintf(stderr,"\t%d unknown\n",nUnknown);
  fflush(NULL);
  return 0;
}
