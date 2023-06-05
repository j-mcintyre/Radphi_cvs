/* map_create_subsystem.c 
   created July 4 1999 by Craig Steffen. */

#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>

void main(int nargs, char *args[]){
  FILE *infile=NULL;
  if(nargs!=3){
    fprintf(stdout,"usage: map_create_subsystem <mapfilename> <new subsystem name>\n");
    exit(1);
  }
  fprintf(stdout,"\n\nCreating new subsystem named %s in map %s\n\n",
	  args[2],args[1]);
  infile=fopen(args[1],"w");
  if(infile==NULL){
    fprintf(stderr,"Problem opening file %s, exiting.\n",args[1]);
    exit(1);
  }
  fprintf(stdout,"Opened map file %s for modification.\n",args[1]);
  if(map_add_sub(infile,args[2],0)){
    fprintf(stderr,"Problem modifying map; exiting.\n");
    exit(1);
  }
  fprintf(stdout,"Creation successful.\n");
}
