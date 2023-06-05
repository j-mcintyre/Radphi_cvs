/*
  $Log: xrdbToMap.c,v $
  Revision 1.4  2004/04/02 17:05:24  radphi
  Makefile - include libdata.a in the link list [rtj]
  xrdb*.c - changed type of main() from void to int to satisfy gcc [rtj]

  Revision 1.3  1997/05/03 15:18:32  radphi
  Changes made by lfcrob@dustbunny
  Added ability to replace files

 * Revision 1.2  1997/02/26  19:16:18  lfcrob
 * Dixed core dump after printing usage
 *
 * Revision 1.1  1997/02/26  19:02:13  lfcrob
 * Initial revision
 *
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <map_manager.h>

static const char rcsid[]="$Id: xrdbToMap.c,v 1.4 2004/04/02 17:05:24 radphi Exp $";

void Usage(void);
int main(int argc, char *argv[]);

void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  xrdbToMap <options> dbFile runNo\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-r\t\tReplace the dbFIle for run runNo.\n");
  fprintf(stderr,"\t-h\t\tPrint this message.\n");
}

int main(int argc, char *argv[])
{
  int iarg;
  char *argptr;
  char *db=NULL;
  char dbFile[1024];
  char *tmp;
  char *fname;
  FILE *fp;
  char *text;
  struct stat buf;
  char map[1024];
  char *subsys="xrdbFiles";
  int size;
  int type;
  int mapSize=1024*1024;
  int runNo;
  int ret;
  char *dir;
  int doReplace=0;

  if((dir=getenv("RADPHI_CONFIG"))==NULL){
    fprintf(stderr,"Environment variable RADPHI_CONFIG not set.\n");
    exit(1);
  }
  sprintf(map,"%s/maps/xrdbFiles.map",dir);

  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'r':
	doReplace=1;
	break;
      case 'h':
	Usage();
	exit(1);
      default:
	fprintf(stderr,"Unknown argument : %s\n",argv[iarg]);
	Usage();
	exit(1);
      }
    }
    else
      break;

  if((argc-iarg) != 2){
    Usage();
    exit(1);
  }

  /* This could be done in a slicker manner..... */

  fname=argv[argc-2];
  strcpy(dbFile,fname);

  /* Now extract the dbName */

  tmp = strtok(dbFile,"/");
  while(tmp != NULL){
    strcpy(dbFile,tmp);
    tmp = strtok(NULL,"/");
  }
  /* Strip off the extension */

  tmp = strtok(dbFile,".");
  strcpy(dbFile,tmp);
  runNo = atoi(argv[argc-1]);
  
  fp = fopen(fname,"r");
  if(!fp){
    fprintf(stderr,"Unable to open %s\n",fname);
    exit(1);
  }
  fprintf(stderr,"Placing %s in %s->%s->%s (run %d)\n",fname,map,subsys,dbFile,runNo);
  fstat(fileno(fp),&buf);
  text = calloc(buf.st_size+1,sizeof(char));
  fread(text,sizeof(char),buf.st_size,fp);
  
  /* Now do the map stuff */

  map_log_mess(0,0);
  if((ret=map_stat_item(map,subsys,dbFile,&size,&type)) != 0){
    switch(ret){
    case MAP_SYSTEM_ERROR_OPEN:
      fprintf(stderr,"%s does not exist.\n",map);
      exit(1);
    case MAP_USER_WARN_NOMATCH:
      fprintf(stderr,"Creating %s->%s->%s\n",map,subsys,dbFile);
      map_add_item(map,subsys,dbFile,mapSize,MAP_CHAR);
      break;
    default:
      exit(1);
    }
  }
  map_log_mess(0,1);
  if(doReplace){
    map_log_mess(MAP_USER_WARN_NOFIND,0);
    ret=map_rem_arr(map,subsys,dbFile,runNo);
    if(ret){
      fprintf(stderr,"%s->%s->%s does not exist for run %d\n",map,subsys,dbFile,runNo);
    }
    map_log_mess(MAP_USER_WARN_NOFIND,1);
  }
  ret=map_put_char(map,subsys,dbFile,buf.st_size+1,text,runNo);
  if(!ret){
    fprintf(stderr,"%s.xrdb placed into %s->%s->%s for run %d\n",fname,map,subsys,dbFile,runNo);
  }
}

