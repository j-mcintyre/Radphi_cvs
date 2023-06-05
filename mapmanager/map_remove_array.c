/* map_create_item.c 
   created Jan 7 2001 by Craig Steffen. */

#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>

int main(int nargs, char *args[]){
  int run_num;
  if(nargs!=5){
    fprintf(stdout,"usage: map_remove_array <mapfilename> <subsystem name> <item> <time>\n");
    exit(1);
  }
  run_num=atoi(args[4]);
  fprintf(stdout,"\n\nRemoving instance %d of item %s in subsystem %s in map %s\n\n",
	  run_num,args[3],args[2],args[1]);
  fprintf(stdout,"Opened map file %s for modification.\n",args[1]);
  if(map_rem_arr(args[1],args[2],args[3],run_num)){
    fprintf(stderr,"Problem modifying map; exiting.\n");
    exit(1);
  }
  fprintf(stdout,"Removal successful.\n");
}
