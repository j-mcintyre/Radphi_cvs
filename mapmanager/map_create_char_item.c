/* map_create_char_item.c 
   created March 8 2001 by Craig Steffen. */

#include <stdio.h>      /* I/O definitions                      */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>

int main(int nargs, char *args[]){
  int length;
  if(nargs!=5){
    fprintf(stdout,"usage: map_create_char_item <mapfilename> <subsystem name> <new item> <length>\n");
    exit(1);
  }
  length=atoi(args[4]);
  fprintf(stdout,"\n\nCreating new item %s in subsystem %s in map %s\n\n",
	  args[3],args[2],args[1]);
  fprintf(stdout,"Opened map file %s for modification.\n",args[1]);
  if(map_add_item(args[1],args[2],args[3],length,MAP_CHAR)){
    fprintf(stderr,"Problem modifying map; exiting.\n");
    exit(1);
  }
  fprintf(stdout,"Creation successful.\n");
}
