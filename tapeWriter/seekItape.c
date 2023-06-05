/*
 * seekItape.c   - seek to a run on the tape
 *
*/

static const char sccsid[] = "@(#)"__FILE__"\t1.7\tCreated 4/4/97 15:08:38, \tcompiled "__DATE__;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <disData.h>
#include <dataIO.h>

#include <mtime.h>

void usage(void)
{
  fprintf(stderr,"Usage: seekItape tapeDevice runNo\n");
  fprintf(stderr,"   runNo can be: '0'-scan the tape and list all encountered runs\n");
  fprintf(stderr,"   runNo can be: >0 -find the run assuming the runs on thetape are sorted (i.e. raw data tapes)\n");
  fprintf(stderr,"   runNo can be: <0 -find the run by seeking the whole tape (i.e. processed data and skimmed tapes)\n");
  exit(1);
}

void main(int argc,char*argv[])
{
  FILE*ifp;
  int runNo = 0;
  int ret;
  int tapeNo = 0;

  if (argc<3) usage();

  ifp = fopen(argv[1],"r");
  if (ifp == NULL)
    {
      fprintf(stderr,"seekItape: Could not open input file %s: %s\n",argv[1],strerror(errno));
      exit(1);
    }

  runNo = strtoul(argv[2],NULL,0);

  ret = data_SeekTape(fileno(ifp),runNo,&tapeNo);

  if (ret) exit(1);

  printf("Tape %d, run %d\n",tapeNo,runNo);

  exit(0);
}

/* end file */
