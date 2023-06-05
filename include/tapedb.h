/*
  $Log: tapedb.h,v $
  Revision 1.1  1997/05/20 09:07:40  radphi
  Initial revision by lfcrob@jlabs2
  Tape database routines

  */

#ifndef TAPEDB_H_INCLUDED
#define TAPEDB_H_INCLUDED

#include <ntypes.h>

typedef struct{
  int runNo;
  int nEvents;
  uint32 bytes;
} tapedbOneRun_t;

typedef struct{
  int tapeNo;
  int nEvents;
  uint32 bytes;
} tapedbOneTape_t;

typedef struct{
  int nRuns;
  tapedbOneRun_t run[1];
} tapedbTape_t;

typedef struct{
  int nTapes;
  tapedbOneTape_t tape[1];
} tapedbRun_t;


int tapedbGetRun(int runNo, tapedbRun_t **run);
int tapedbGetTape(int tapeNo, tapedbTape_t **tape);
int tapedbPut(int tapeNo,int runNo,int nEvents, uint32 bytes);
int tapedbInit(void);
#endif
