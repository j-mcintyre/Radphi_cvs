/*
  $Log: dumpDirToTapes.c,v $
  Revision 1.1  1997/07/31 13:44:18  radphi
  Initial revision by lfcrob@jlabs1
  Dump a directory to one or more tapes

  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <signal.h>

#include <tapeIO.h>
#include <tapeData.h>
#include <dataIO.h>
#include <disData.h>
#include <cefIO.h>
#include <tapedb.h>

#define BUFSIZE 100000

void sigHandle(int signal);


typedef struct{
  int nEvents;
  int kBytes;
  int nPartial;
  int nComplete;
} tapeSession_t;
typedef struct{
  char name[1024];
  int runNo;
} file_t;


void main(int argc, char *argv[]);
int sortFiles(const void *one,const void *two);

static int running=1;

void main(int argc, char *argv[])
{
  int iarg;
  char *argptr;
  char *directory=NULL;
  char *outputs[10];
  FILE *fpOut[10];
  tapeSession_t tapeInfo[10];
  tapeSession_t tapeSessionInfo[10];
  int onTape[10];
  int tapeNo[10];
  int needToWrite;
  int nOutputs=0;
  int index,j;
  struct dirent *entry;

  DIR *dir;
  char fullPath[1024];
  struct stat fStatus;
  itape_header_t *event;
  int *cefBuf;
  int cefHandle;
  int nProcessed=0;
  int nEvents;
  int nBytes;
  int nBytesProcessed=0;
  int ret;
  int runNo;
  int i=0;

  file_t files[1024];
  int nFiles=0;


  tapedbTape_t *tapedb[10];
  tapedbRun_t *rundb;

  event=malloc(BUFSIZE);
  cefBuf = malloc(BUFSIZE);
  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'd':
	directory = argv[++iarg];
	break;
      }
    }
    else{
      outputs[nOutputs++] = argv[iarg];
    }

  if(!directory)
    directory = ".";
  fprintf(stderr,"Dumping directory : %s\n",directory);
  
  for(index=0;index<nOutputs;index++){
    fprintf(stderr,"\t%s\n",outputs[index]);    
  }



  /* Open the tapes, check headers */

  
  for(index=0;index<nOutputs;index++){
    tapeHeader_t *tHeader;
    int first=1;

    fpOut[index] = fopen(outputs[index],"r");
    if(!fpOut[index]){
      fprintf(stderr,"Could not open %s for reading!\n",outputs[index]);
      exit(0);
    }
    tape_rewind(fileno(fpOut[index]));
    if(data_read(fileno(fpOut[index]),event,BUFSIZE) != DATAIO_OK){
      fprintf(stderr,"Failed to read tape header from %s\n",outputs[index]);
      exit(0);
    }
    if(event->type != TYPE_TAPEHEADER){
      fprintf(stderr,"No tape header on %s\n",outputs[index]);
      exit(0);
    }
    tHeader= (void *) event;
    tapeNo[index] = tHeader->tapeSeqNumber;
    fprintf(stderr,"Writing tape %d in %s\n",tHeader->tapeSeqNumber,outputs[index]);    
    
    memset(&tapeInfo[index],0,sizeof(tapeSession_t));
    memset(&tapeSessionInfo[index],0,sizeof(tapeSession_t));
    if(!tapedbGetTape(tHeader->tapeSeqNumber,&tapedb[index])){
      for(runNo=0;runNo<tapedb[index]->nRuns;runNo++){
	if(tapedb[index]->run[runNo].nEvents == -1)
	  tapeInfo[index].nPartial++;
	else{
	  tapeInfo[index].nComplete++;
	  tapeInfo[index].nEvents += tapedb[index]->run[runNo].nEvents;
	  tapeInfo[index].kBytes += tapedb[index]->run[runNo].bytes/1024;
	}
      }
      fprintf(stderr,"Tape %d already has %d kBytes on it in %d runs  (%d complete, %d partial)\n",tapeNo[index],
	      tapeInfo[index].kBytes,
	      tapeInfo[index].nPartial,
	      tapeInfo[index].nComplete,
	      tapeInfo[index].nPartial-tapeInfo[index].nComplete);
    }
    else{
      tapedb[index]=NULL;
    }


    /* Read the tape database for this tape */

    ret = 1;
    while(ret>=0){
      if(first){
	tape_flush(fileno(fpOut[index]));
	fprintf(stderr,"Seek to EOD");
	fflush(stderr);
	first=0;
      }
      fprintf(stderr,".");
      fflush(stderr);
      ret=tape_findFM(fileno(fpOut[index]),10);
    }
    fprintf(stderr,"at EOD\n");

    tape_flush(fileno(fpOut[index]));
    fclose(fpOut[index]);
    fpOut[index] =  fopen(outputs[index],"w");
  }

  /* Open the directory */
  
  dir = opendir(directory);
  if(!dir){
    fprintf(stderr,"Failed to open %s\n",directory);
    exit(1);
  }

  /* Read and Sort the directory */

  while(((entry = readdir(dir)) != NULL)){
    sprintf(fullPath,"%s/%s",directory,entry->d_name);
    stat(fullPath,&fStatus);
    if(fStatus.st_mode & S_IFREG){
      if((strstr(entry->d_name,".evt") != NULL) && (strstr(entry->d_name,"run") != NULL)){
	strtok(entry->d_name,"_");
	runNo = atoi(strtok(NULL,"_"));
	strcpy(files[nFiles].name,fullPath);
	files[nFiles].runNo = runNo;
	nFiles++;
      }
    }
  }
  qsort(files,nFiles,sizeof(file_t),sortFiles);
  

  signal(SIGHUP,sigHandle);
  signal(SIGINT,sigHandle);
  i=0;
  while(running && (i < nFiles)){
    sprintf(fullPath,"%s",files[i].name);
    stat(fullPath,&fStatus);
    runNo = files[i].runNo; 
    if(fStatus.st_mode & S_IFREG){
	
      /* Check to make sure this run is not on the tape already */
    
      needToWrite=0;
      for(index=0;index<nOutputs;index++){
	onTape[index]=0;
	if(tapedb[index]){
	  for(j=0;j<tapedb[index]->nRuns;j++){
	    if(tapedb[index]->run[j].runNo == runNo){
	      onTape[index]=1;
	      fprintf(stderr,"Run %d already exists on %s\n",runNo,outputs[index]);
	      break;
	    }
	  }
	  if(!onTape[index]){
	    needToWrite=1;
	  }
	}
	else{
	  needToWrite=1;
	}
      }
      
      if(!needToWrite){
	fprintf(stderr,"Run %d from %s already exists for all output streams, skipping\n",runNo,fullPath);
      }
      
      if(fStatus.st_size == 0)
	fprintf(stderr,"Run %d from %s skipped (empty file)\n",runNo,fullPath);
      
      /* Open the CEF file */
      
      if((fStatus.st_size != 0) && (needToWrite)){
	fprintf(stderr,"Writing %d kBytes from %s\n",fStatus.st_size,fullPath);
	if(evOpen(fullPath,"r",&cefHandle) == S_SUCCESS){
	  
	  nEvents=0;
	  nBytes=0;
	  for(index=0;index<nOutputs;index++){
	    if(!onTape[index]){
	      tapedbPut(tapeNo[index],runNo,-1,0);
	      tapeInfo[index].nPartial++;
	      tapeSessionInfo[index].nPartial++;
	    }		
	  }
	  
	  while(running && (evRead(cefHandle,cefBuf,BUFSIZE) == S_SUCCESS)){
	    if(!wrapEvent(cefBuf,event,BUFSIZE)){
	      data_clean(event);
	      nEvents++;
	      nBytes += event->length;
	      nBytesProcessed += event->length;
	      for(index=0;index<nOutputs;index++)
		if(!onTape[index]){
		  ret = data_write(fileno(fpOut[index]),event);
		  if(ret != DATAIO_OK){
		    fprintf(stderr,"Error writing to tape %s, error = %d\n",outputs[index],ret);
		    running=0;
		  }
		  else{
		    tapeInfo[index].nEvents++;
		    tapeSessionInfo[index].nEvents++;
		  }
		}
	      if((++nProcessed%10) == 0){
		fprintf(stderr,"%d %d\r",nProcessed,nEvents);
		fflush(stderr);
	      }
	    }
	    else{
	      fprintf(stderr,"Failed to wrap event!\n");
	    }
	  }
	  for(index=0;index<nOutputs;index++){
	    if(!onTape[index]){
	      tapeInfo[index].nComplete++;
	      tapeSessionInfo[index].nComplete++;
	      data_flush(fileno(fpOut[index]));
	      tape_writeFM(fileno(fpOut[index]));
	      tapeInfo[index].kBytes += (nBytes/1024);
	      tapeSessionInfo[index].kBytes += (nBytes/1024);
	      if(tapedbPut(tapeNo[index],runNo,nEvents,nBytes)){
		fprintf(stderr,"Fatal error! Connot update database\n");
		exit(1);
	      }
	    }
	  }
	  evClose(cefHandle);
	}
	else{
	  fprintf(stderr,"Could not open %s for reading.\n",fullPath);
	}	
      }
    }
    i++;
  }
  fprintf(stderr,"\n\n*************************************************\n");
  fprintf(stderr,"Processed %d events (%d kBytes)\n",nProcessed,nBytesProcessed);
  for(index=0;index<nOutputs;index++){
    fprintf(stderr,"\n");
    fprintf(stderr,"Device: %s\n",outputs[index]);
    fprintf(stderr,"Tape Number: %d\n",tapeNo[index]);
    fprintf(stderr,"\tSession: %7d kBytes %4d runs %7d events\n",tapeSessionInfo[index].kBytes,tapeSessionInfo[index].nPartial,tapeSessionInfo[index].nEvents);
    fprintf(stderr,"\tTotal  : %7d kBytes %4d runs %7d events\n",tapeInfo[index].kBytes,tapeInfo[index].nPartial,tapeInfo[index].nEvents);
  }
  
}

void sigHandle(int signal)
{
  fprintf(stderr,"\nInterrupt!\n");
  running=0;
}
int sortFiles(const void *a,const  void *b)
{
  file_t *one = (file_t *)a;
  file_t *two = (file_t *)b;
  if(one->runNo < two->runNo)
    return(-1);
  if(one->runNo > two->runNo)
    return(1);
  return(0);
      
}
