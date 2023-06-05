/*
 * feedFromFile.c  - feed the dispatcher from an itape file
*/

static const char sccsid[] = "@(#)"__FILE__"\t1.30\tCreated 4/29/96 17:18:16, \tcompiled "__DATE__;

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <mtime.h>
#include <math.h>
#include <limits.h>
#ifndef HPUX
#include <sys/select.h>
#endif

#include <disData.h>
#include <tapeData.h>
#include <dataIO.h>
#include <tapeIO.h>
#include <disIO.h>
#include <wrapper.h>
#include <msql.h>

/* DD stuff */

#include <semaphore.h>
#ifndef _LANGUAGE_C
#define _LANGUAGE_C
#endif
#include <dd_user.h>


#define BUFSIZE 100000

void usage(void)
{
  fprintf(stderr,"Usage: feedFromFile switches... files...\n");
  fprintf(stderr,"Switches are:\n");
  fprintf(stderr,"\t-D<hostname:port>\tconnect to the dispatcher at host:port\n");
  fprintf(stderr,"\t-H<port>\tuse this dispatcher TCP port\n");
  fprintf(stderr,"\t-M<nevents>\tmaximum number of events to send\n");
  fprintf(stderr,"\t-p\t\tsend processed data\n");
  fprintf(stderr,"\t-t<tapeNo>\tuse <tapeNo> as tape number if the input data has no tape labels\n");
  fprintf(stderr,"\t-r<runNo>\tseek the tape to <runNo> and process only this run, give a negative run number to disable 'smart' seeking.\n");
  fprintf(stderr,"\t-l<latch>\tset the latch mask to check trigger=1 events against\n");
  fprintf(stderr,"\t-S<spillNo>\tsend data starting from this spill\n");
  fprintf(stderr,"\t-k<n> \t\tlimit data rate to n kb/sec\n");
  fprintf(stderr,"\t-e<n> \t\tlimit data rate to n events/sec\n");
  fprintf(stderr,"\t-R\t\trewind tape when finished\n");
  fprintf(stderr,"\t-E\t\teject tape when finished\n");
  fprintf(stderr,"\t-x<number>\tenable read error recovery strategy <number>. Default: 0\n");
  fprintf(stderr,"\t\t\t-x0 - no recovery, exit on errors\n");
  fprintf(stderr,"\t\t\t-x1 - flush the input buffer and retry the read request. Exit after 5 consecutive errors.\n");
  fprintf(stderr,"\t-h \t\tprint this message\n");

  exit(1);
}

static int FirstSpill = 0;

static int crunNo = 0;

static uint32 outputLength = 0;

static int eventsThisRun = 0;
static int timeThisRun = 0;
static int eventsThisRunSent = 0;
static int eventsThisTape = 0;
static int eventsThisTapeSent = 0;

static int tapeNo = 0;

static int haveBeginTape = 0;

static int exitFlag = 0;

static int requestedRaw = 0;

int waitForCommand(const char* cmd)
{
/*  fprintf(stderr,"DisFeed: Waiting for command: %s\n",cmd);*/

  while (1)
    {
      int ret;

      fd_set rset;
      int maxfd = 0;

      FD_ZERO(&rset);
      FD_SET(disIO_socket,&rset); if (disIO_socket > maxfd) maxfd = disIO_socket;

      ret = select(maxfd + 1, &rset, NULL, NULL, NULL);
      
      if (ret < 0)
	{
	  fprintf(stderr,"DisFeed: Error: select() returned %d, errno: %d (%s)\n",ret,errno,strerror(errno));
	  return -1;
	}

      if (FD_ISSET(disIO_socket,&rset))
	{
	  /* ok, we have something to read */

	  static void *msg = NULL;
	  static int msglen = 0;
	  char *cmd0;
	  char *word;
	  int maybeUnexpected = 0;

	  if (msg)
	    free(msg);

	  msg = NULL;

	  ret = disIO_readRAW_alloc(&msg,&msglen,0);

	  if (ret == DISIO_EOF)
	    {
	      fprintf(stderr,"DisFeed: Error: Lost connection to the Dispatcher.\n");
	      return -1;
	    }

	  if (ret != DISIO_COMMAND)
	    {
	      fprintf(stderr,"DisFeed: Error: Unexpected message from Dispatcher: %d, expected a DISIO_COMMAND (%d)\n",
		      ret, DISIO_COMMAND);

	      continue;
	    }

/*	  fprintf(stderr,"DisFeed: Command from Dispatcher: %s\n",msg);*/

	  /* ok, we got a command. Now parse it */

	  word = strtok(msg,":");
	  cmd0 = word;

	  if (strcmp(word,"NOP") == 0)
	    {
	      /* do nothing */
	    }
	  else if (strcmp(word,"PING") == 0)
	    {
	      disIO_command("PING-REPLY");
	    }
	  else if (strcmp(word,"REQUEST_DATA") == 0)
	    {
	      /* fprintf(stderr,"DisFeed: Command from Dispatcher: %s\n",word); */

	      word = strtok(NULL,":");

	      if (word)
		requestedRaw += strtoul(word,NULL,0);
	      else
		requestedRaw += 1;
	    }
	  else if (strcmp(word,"QUIT") == 0)
	    {
	      fprintf(stderr,"DisFeed: Command from Dispatcher: %s\n",word);

	      exitFlag = 1;
	    }
	  else
	    {
	      maybeUnexpected = 1;
	    }

	  if (strcmp(cmd0,cmd) == 0)
	    {
	      maybeUnexpected = 0;

	      /* fprintf(stderr,"DisFeed:  Got what we waited for: %s\n",word); */

	      return 0;
	    }

	  if (maybeUnexpected)
	    fprintf(stderr,"DisFeed: Unexpected command from Dispatcher: '%s' was ignored.\n",msg);
	}

    }
}

int sendEndtape(void)
{
  char reply[256];
  int nread;
  int ret;
  unsigned long currentTime = mtime();
  unsigned long runTime = currentTime - timeThisRun;

  printf("DisFeed: Sending ENDTAPE...\n");
  disIO_command("ENDTAPE");
  
  printf("DisFeed: Endrun %d: Total events: %d/%d, Total time: %d seconds (%6.2lf hours), EventsPerSecond: %6.2lf\n",
	 crunNo,eventsThisRunSent,eventsThisRun,(int)(runTime/1000.0),(double)((runTime/1000.0)/(double)3600),(double)(eventsThisRun/(double)(runTime/1000.0)));
  printf("  Total events this tape: %d/%d\n",eventsThisTapeSent,eventsThisTape);

  printf("DisFeed: Waiting for ALLFINISHED...\n");

  ret = waitForCommand("ALLFINISHED");

  if (ret != 0)
    {
      fprintf(stderr,"DisFeed: Error: Have problems talking to the Dispatcher, exiting...\n");
      exit(1);
    }
  
  haveBeginTape = 0;
  
  return 0;
}

double limitKbps=0;
double limitEvps=0;

void main(int argc,char*argv[])
{
  int i;
  char*arg;
  char myName[128];
  char *host = getenv("RADPHI_DISPATCHER");
  int sendProcessed = 0;
  int eofCount = 0;
  int runNo = 0;
  uint32 latch = 0xffffffff;
  int port =10099;
  int MaxEvents = 0;
  int nevents = 0;
  int wantTapeNo = 0;
  int doRewind = 0;
  int doEject = 0;
  int recover = 0;
  int maxReadRetries = 5;
  int readRetry;
  itape_header_t *ih=NULL;

  /* DD stuff */

  fifo_mode fmode;
  struct fifo_entry fev;
  int go_on=1;
  int status;
  int ctl[4]={-1,-1,-1,-1};

  fmode.mode = FMODE_ONREQ;
  fmode.wait = DD_WAIT_SLEEP;
  fmode.prescale=1;
  fmode.suser=0;
  fmode.p2ctl = ctl;


  setbuf(stdout,NULL);
  setbuf(stderr,NULL);

  printf("DisFeed: Command line arguments:\n");
  for (i=0; i<argc; i++)
    {
      printf("DisFeed: argv[%d]: '%s'\n",i,argv[i]);
    }
  printf("DisFeed: End of Command line arguments\n");

  for (i=1; i<argc; i++)
    {
      if ((*argv[i])=='-')
	{
	  arg = &argv[i][2];
	  switch (argv[i][1])
	    {
	    case 'k':
	      limitKbps = atof(arg);
	      if(limitKbps < 0)
		limitKbps=0;
	      break;
	    case 'e':
	      limitEvps = atof(arg);
	      if(limitEvps < 0)
		limitEvps=0;
	      break;
	    case 'D': 
	      host = arg; 
	      break;
	    case 'H': 
	      port = atoi(arg); 
	      break;
	    case 'r': 
	      runNo = atoi(arg); 
	      break;
            case 'l':
              latch = strtoul(arg,NULL,0);
              break;
	    case 't': 
	      wantTapeNo = (int)strtoul(arg,NULL,0); 
	      break;
	    case 'p': 
	      sendProcessed = 1; 
	      break;
	    case 'S': 
	      FirstSpill = atoi(arg);
	      break;
	    case 'M': 
	      MaxEvents = (int)strtoul(arg,NULL,0);
	      break;
	    case 'R': 
	      doRewind = 1;
	      break;
	    case 'E': 
	      doEject = 1;
	      break;
	    case 'x':
	      recover = (int)strtoul(arg,NULL,0);
	      break;
	    case 'h': 
	      usage(); 
	      break;
	    default:
	      fprintf(stderr,"DisFeed: Unknown switch %s\n",argv[i]);
	      usage();
	      break;
	    }
	}
    }

  {
    int ret;
    ret = disIO_connect(host,port);
    if (ret)
      {
	fprintf(stderr,"DisFeed: Cannot connect to the Dispatcher on %s, exiting...\n",host);
	exit(1);
      }
  }

  exitFlag = 0;

  /* Open DD connection */

  ddu_init("Dispatcher",fmode);
  ddu_start();

  fprintf(stderr,"Connected to the DD ring\n");


  /* Set default run number and time */

  {
    int sock;
    m_result *result;
    m_row row;
    char *query="SELECT runNumber FROM sessions";
    time_t secs;
    int runNo;
    /* Get from msql */

    fprintf(stderr,"Getting run number from mSQL\n");
    sock = msqlConnect("urs1.cebaf.gov");
    msqlSelectDB(sock,"radphi_exp");
    msqlQuery(sock,query);
    result = msqlStoreResult();
    row = msqlFetchRow(result);
    runNo = atoi(*row);    
    fprintf(stderr,"Run number %d\n",runNo);
    msqlFreeResult(result);    
    secs = time(NULL);
    wrapperSetDefaults(WRAP_DEF_RUN,runNo);
    wrapperSetDefaults(WRAP_DEF_TIME,secs);
  }


  while(exitFlag==0 && ddu_attached()){
    void *buf;
    int ret;
    
    if(ih == NULL){
      ih = malloc(BUFSIZE);
    }

    ret = ddu_get_fev(&fev);
    if(ret!=0){
      fprintf(stderr,"ddu_get_fev() failed, error code %d\n",ret);
      exitFlag=1;
      continue;
    }
    else{

      /* Wrap the event */

      if(wrapEvent(fev.p2da,ih,BUFSIZE)){
	fprintf(stderr,"Failed to wrap event!\n");
	exitFlag=1;
	continue;
      }
      if((ret = ddu_put_fev(fev)) != 0){
        fprintf(stderr,"ddu_put_fev() returned %d\n",status);
        exitFlag=1;
	continue;
      }
    }
	      
    {   /*Main*/
      static int cspillNo = 0;
      
      if (ih->runNo != crunNo){
	char msg[256];
	
	printf("DisFeed: Run: %d\n",ih->runNo);
	
	if (crunNo != 0)
	  {
	    sendEndtape();
	  }
	
	if (runNo > 0)
	  if (ih->runNo != runNo)
	    {
	      fprintf(stderr,"DisFeed: Run number mismatch: expected %d, read %d, exiting...\n",
		      runNo,ih->runNo);
	      exitFlag = 1;
	      break;
	    }
	
	if (tapeNo == 0) tapeNo = wantTapeNo;
	
	sprintf(msg,"BEGINTAPE:%d:%d",ih->runNo,tapeNo);
	
	printf("DisFeed: Sending BEGINTAPE...\n");
	disIO_command(msg);
	
	haveBeginTape = 1;
	
	crunNo = ih->runNo;
	cspillNo = 0;
	
	eventsThisTape += eventsThisRun;
	eventsThisTapeSent += eventsThisRunSent;
	eventsThisRun = 0;
	eventsThisRunSent = 0;
	timeThisRun = (int)mtime();
	
      } /* endif (begin a new run) */
      
      if (ih->eventNo%100 == 0){
	static int stotal = 0;
	static unsigned long stime = 0;
	static int stotalSent = 0;
	static int count = 0;
	double epsSent = 0;
	double eps = 0;
	double kbps = 0;
	unsigned long currentTime = mtime();
	static unsigned long cntLength = 0;
	
	if (stime == 0) stime = timeThisRun;
	
	if (eventsThisRun - stotal > 0)
	  {
	    eps = (double)(eventsThisRun - stotal)/(double)((currentTime - stime)/1000.0);
	    epsSent = (double)(eventsThisRunSent - stotalSent)/(double)((currentTime - stime)/1000.0);
	    kbps = (double)(outputLength-cntLength)/(double)((currentTime-stime)/1000.0)/1024;
	    
	    if(limitKbps != 0)
	      if(finite(kbps)){
		unsigned long naptime;
		unsigned int usleepnaptime; /* in microsec */
		if(kbps > limitKbps){
		  naptime = (unsigned long)((outputLength-cntLength)/(1024.0*limitKbps/1000.0)) + (unsigned long)stime;  /* millisec to nap + stime */
		  usleepnaptime = (int)((naptime-currentTime)*1000);
#ifdef HPUX
		  while(1){
		    if(mtime() > naptime)
		      break;
		  }
#else
		  usleep(usleepnaptime);
#endif
		}
	      }
	    if(limitEvps != 0)
	      if(finite(kbps)){
		unsigned long naptime;
		unsigned int usleepnaptime; /* in microsec */
		if(kbps > limitEvps){
		  naptime = (unsigned long)((eventsThisRunSent - stotalSent)/(limitEvps/1000.0)) + stime; /* millisec to nap + stime */
		  usleepnaptime = (int)((naptime-currentTime)*1000);
#ifdef HPUX
		  while(1){
		    if(mtime() > naptime)
		      break;
		  }
#else
		  usleep(usleepnaptime);
#endif
		}
	      }
	    currentTime = mtime();
	    eps = (double)(eventsThisRun - stotal)/(double)((currentTime - stime)/1000.0);
	    epsSent = (double)(eventsThisRunSent - stotalSent)/(double)((currentTime - stime)/1000.0);
	    kbps = (double)(outputLength-cntLength)/(double)((currentTime-stime)/1000.0)/1024;
	  }
	
	count++;
	if (count>10)
	  {
	    time_t timelog;
	    time(&timelog);
	    printf("Time stamp: %s",ctime(&timelog));
	    count = 0;
	  }
	
	printf("DisFeed: Run: %d, Spill: %3d, TotalEvents: %7d/%7d, Events/Second: %6.2lf/%6.2lf, Kbytes/sec: %10.3lf\n",ih->runNo,ih->spillNo,eventsThisRunSent,eventsThisRun,epsSent,eps,kbps);
	cntLength = outputLength;
	
	stotal = eventsThisRun;
	stotalSent = eventsThisRunSent;
	stime = currentTime;
	
	cspillNo = ih->spillNo;
	
      } /* endif (begin a new spill) */
    }
    
    if (ih->spillNo > FirstSpill) {
      FirstSpill = 0;
    }
    
    if (!FirstSpill) 
      {
	eventsThisRun ++;
	nevents ++;
	outputLength += ih->length + 4; /* extra 4 bytes of CRC per event */
	if (!(ih->trigger==1)||!(ih->latch)||(latch&ih->latch))
	  {
	    if(!sendProcessed){
	      
	      if (requestedRaw <= 0)
		{
		  waitForCommand("REQUEST_DATA");
		}
	      
	      requestedRaw --;
	    }
	    
	    if (sendProcessed)
	      ret = disIO_writeCOOKED(ih,ih->length + 4);
	    else
	      ret = disIO_writeRAW(ih,ih->length + 4);
	    
	    switch (ret)
	      {
	      case DISIO_OK:
		eventsThisRunSent++;
		break;
	      case DISIO_EOF:
		fprintf(stderr,"DisFeed: Error: Dispatcher disconnected, exiting...\n");
		exitFlag = 1;
		break;
	      default:
		fprintf(stderr,"DisFeed: Error: disIO_write() return code %d, exiting...\n",ret);
		exitFlag = 1;
		break;
	      }
	  }
      }


    if (MaxEvents) 
      if (nevents > MaxEvents) 
	exitFlag = 1;

    
  } /* end master While */
  

  if (haveBeginTape){
    sendEndtape();
  }

  printf("DisFeed: Disconnecting from the Dispatcher...\n");
  disIO_disconnect();
  
  ddu_stop();
  ddu_close();


  printf("DisFeed: Finished.\n");
  exit(0);
  
  
}

/* end of file */
