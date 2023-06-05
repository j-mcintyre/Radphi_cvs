/*
  $Log: xrdbFromMap.c,v $
  Revision 1.2  2004/04/02 17:05:24  radphi
  Makefile - include libdata.a in the link list [rtj]
  xrdb*.c - changed type of main() from void to int to satisfy gcc [rtj]

  Revision 1.1  1997/02/26 19:02:34  lfcrob
  Initial revision

  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <map_manager.h>

static const char rcsid[]="$Id";


int main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
  char *dir;
  char map[1024];
  char data[1024*1024];
  char *itemName;
  int runNo;
  int olength;
  int actualRun;

  if((dir=getenv("RADPHI_CONFIG"))==NULL){
    fprintf(stderr,"Environment variable RADPHI_CONFIG not set.\n");
    exit(1);
  }
  sprintf(map,"%s/maps/xrdbFiles.map",dir);
  if(argc != 3){
    fprintf(stderr,"Usage : xrdbFromMap itemName runNo\n");
    exit(1);
  }
  itemName = argv[1];
  runNo = atoi(argv[2]);
  fprintf(stderr,"Requested:\n");
  fprintf(stderr,"%s->xrdbFiles->%s->%d\n",map,itemName,runNo);
  map_get_char(map,"xrdbFiles",itemName,sizeof(data),data,runNo,&actualRun,&olength);
  fprintf(stderr,"Returned:\n");
  fprintf(stderr,"%s->xrdbFiles->%s->%d\n",map,itemName,actualRun);
  printf("%s\n",data);
  
}
