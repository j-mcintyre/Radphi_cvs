/*
  $Log: tapedb.c,v $
  Revision 1.1  1997/05/20 09:08:05  radphi
  Initial revision by lfcrob@jlabs2
  Tape database routines

  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tapedb.h>

static char *tapedbDir=NULL;


int tapedbInit(void)
{
  if(!tapedbDir){
    if((tapedbDir = getenv("RADPHI_TAPEDB")) == NULL){
      fprintf(stderr,"tapedbInit: RADPHI_TAPEDB not defined\n");
      return(1);
    }    
  }
  return(0);
}



int tapedbPut(int tapeNo,int runNo,int nEvents, uint32 bytes)
{
  FILE *fpRun;
  char runFile[1024];
  FILE *fpTape;
  char tapeFile[1024];

  if(!tapedbDir){
    if(tapedbInit())
      return(1);
  }
  sprintf(runFile,"%s/run_%d.tdb",tapedbDir,runNo);
  if((fpRun = fopen(runFile,"a")) == NULL)
    return(1);
  sprintf(tapeFile,"%s/tape_%d.tdb",tapedbDir,tapeNo);
  if((fpTape = fopen(tapeFile,"a")) == NULL)
    return(1);
   
  /* Write the date */
  fprintf(fpRun,"%d %d %d\n",tapeNo,nEvents,bytes);
  fprintf(fpTape,"%d %d %d\n",runNo,nEvents,bytes);
  fclose(fpRun);
  fclose(fpTape);
  return(0);
}


int tapedbGetTape(int tapeNo, tapedbTape_t **tape)
{
  int nRuns=0;
  int runNo;
  int nEvents;
  int bytes;
  FILE *fpTape;
  char tapeFile[1024];

  if(!tapedbDir){
    if(tapedbInit())
      return(1);
  }
  sprintf(tapeFile,"%s/tape_%d.tdb",tapedbDir,tapeNo);
  if((fpTape = fopen(tapeFile,"r")) == NULL)
    return(1);
  while(fscanf(fpTape,"%d %d %d",&runNo,&nEvents,&bytes) == 3){
    nRuns++;
  }

  

  *tape = malloc(sizeof(tapedbTape_t) + (nRuns-1)*sizeof(tapedbOneRun_t));
  rewind(fpTape);
  nRuns=0;
  while(fscanf(fpTape,"%d %d %d",&runNo,&nEvents,&bytes) == 3){
    (*tape)->run[nRuns].runNo = runNo;
    (*tape)->run[nRuns].nEvents = nEvents;
    (*tape)->run[nRuns].bytes = bytes;
    nRuns++;
  }
  (*tape)->nRuns = nRuns;
  return(0);
}

int tapedbGetRun(int runNo, tapedbRun_t **run)
{
  int nTapes=0;
  int tapeNo;
  int nEvents;
  int bytes;
  FILE *fpRun;
  char runFile[1024];

  if(!tapedbDir){
    if(tapedbInit())
      return(1);
  }
  sprintf(runFile,"%s/run_%d.tdb",tapedbDir,runNo);
  if((fpRun = fopen(runFile,"r")) == NULL)
    return(1);
  while(fscanf(fpRun,"%d %d %d",&tapeNo,&nEvents,&bytes) == 3){
    nTapes++;
  }

  *run = malloc(sizeof(tapedbRun_t) + (nTapes-1)*sizeof(tapedbOneTape_t));
  rewind(fpRun);
  nTapes=0;
  while(fscanf(fpRun,"%d %d %d",&tapeNo,&nEvents,&bytes) == 3){
    (*run)->tape[nTapes].tapeNo = tapeNo;
    (*run)->tape[nTapes].nEvents = nEvents;
    (*run)->tape[nTapes].bytes = bytes;
    nTapes++;
  }
  (*run)->nTapes = nTapes;
  return(0);
}

