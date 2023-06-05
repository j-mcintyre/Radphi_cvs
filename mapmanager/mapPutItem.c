/*
  $Log: mapPutItem.c,v $
  Revision 1.2  2000/01/18 04:15:26  radphi
  Changes made by jonesrt@golem
  major upgrade to mapamanger package to support little-endian processors;
  data on disk remains in big-endian format as before, and is translated on
  input and output for little-endian machines (such as the Intel x86).
  -rtj

  Revision 1.1  1997/04/08 21:20:33  radphi
  Initial revision

  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <map_manager.h>

static char rcsid[] = "$Id: mapPutItem.c,v 1.2 2000/01/18 04:15:26 radphi Exp $";


int main(int argc, char *argv[]);
void Usage(void);



void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  mapPutItem <options> map subsystem item time\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-r\t\t\tReplace item for given time\n");
  fprintf(stderr,"\t-c <length> <type>\tCreate item of array size <length> and\n\t\t\t\ttype <type> if it does not exist\n");
  fprintf(stderr,"\t-h\t\t\tPrint this message\n");
}


int main(int argc, char *argv[])

{
  int iarg;
  char *argptr;
  int *intArray;
  float *floatArray;
  char  *charArray;

  char *map;
  char *subsystem;
  char *item;
  char *typeChar;
  int type;
  int length;
  int time;

  int doReplace=0;
  int doCreate=0;

  int ret;
  int itemLength;
  int itemType;
  int itemTime;

  int reading=1;
  int nRead=0;


  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'r':
	doReplace=1;
	break;
      case 'c':
	doCreate=1;
	length = atoi(argv[++iarg]);	
	if(length<1){
	  fprintf(stderr,"Length must be >= 1\n");
	  exit(1);
	}
	typeChar = argv[++iarg];
	if(!strcmp(typeChar,"int")){
	  type = MAP_INT;
	}
	else if (!strcmp(typeChar,"float")) {
	  type = MAP_FLOAT;
	}
	else if (!strcmp(typeChar,"char")) {
	  type = MAP_CHAR;
	}
	else{
	  fprintf(stderr,"Unknown type : %s\n",typeChar);
	  Usage();
	  exit(1);      
	}	
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
    else{
      break;
    }
  if(argc-iarg != 4){
    Usage();
    exit(1);
  }
  map = argv[iarg++];
  subsystem = argv[iarg++];
  item = argv[iarg++];
  time = atoi(argv[iarg]);

  map_log_mess(MAP_USER_WARN_NOMATCH,0);
  ret = map_stat_item(map,subsystem,item,&itemLength,&itemType);  
  if(doCreate){
    if(ret == MAP_USER_WARN_NOMATCH){
      /* Ok, we can create it */
      map_add_item(map,subsystem,item,length,type);
      printf("Created %s->%s in map %s\n",subsystem,item,map);
    }
    else{
      /* Already exists */
      fprintf(stderr,"%s->%s already exists in map %s\n",subsystem,item,map);
      exit(1);
    }
  }
  else{
    if(ret == MAP_USER_WARN_NOMATCH){
      fprintf(stderr,"%s->%s does not exist in map %s\n",subsystem,item,map);
      exit(1);
    }
    else{      
      length = itemLength;
      type = itemType;
    }
  }
  map_log_mess(MAP_USER_WARN_NOMATCH,1);

  /* Ok, the subsystem exists - now we must allocate space for data */

  switch(type){
  case MAP_INT:
    intArray = calloc(length,sizeof(int));
    break;
  case MAP_FLOAT:
    floatArray = calloc(length,sizeof(float));
    break;
  case MAP_CHAR:
    charArray = calloc(length,sizeof(char));
    break;
  }

  /* Now see if data exists for that time already */

  map_log_mess(MAP_USER_WARN_NOFIND,0);
  switch(type){
  case MAP_INT:
    ret = map_get_int(map,subsystem,item,length,intArray,time,&itemTime);
    break;
  case MAP_FLOAT:
    ret = map_get_float(map,subsystem,item,length,floatArray,time,&itemTime);
    break;
  case MAP_CHAR:
    exit(1);
    break;
  }
  if(ret == MAP_USER_WARN_NOFIND){
    itemTime = time+1;
  }
  else
    if((time == itemTime) && !doReplace){
      fprintf(stderr,"Data already exists for time %d\n",time);
      exit(1);
    }
  map_log_mess(MAP_USER_WARN_NOFIND,1);


  printf("Enter up to %d white space separated values.\nArray will be filled out with zeros:\n",length);
  while(reading){
    char str[1024];
    ret = scanf("%s",str);
    switch(ret){
    case 1:
      switch(type){
      case MAP_INT:
	intArray[nRead] = strtol(str,NULL,0);
	break;
      case MAP_FLOAT:
	floatArray[nRead] = atof(str);
	break;
      case MAP_CHAR:
	exit(1);
	break;
      }
      nRead++;
      break;
    case 0:
      break;
    case EOF:
      reading=0;
      break;
    default:
      exit(1);
    }
    if(nRead == length)
      reading=0;
  }
  printf("Read %d values\n",nRead);
  if((time == itemTime) && doReplace){
    map_rem_arr(map,subsystem,item,time);
  }
  switch(type){
  case MAP_INT:
    map_put_int(map,subsystem,item,length,intArray,time);
    break;
  case MAP_FLOAT:
    map_put_float(map,subsystem,item,length,floatArray,time);
    break;
  case MAP_CHAR:
    map_put_char(map,subsystem,item,length,charArray,time);
    break;
  }
}
