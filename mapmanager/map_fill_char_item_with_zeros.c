#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <map_manager.h>

static char rcsid[] = "$Id: map_fill_char_item_with_zeros.c,v 1.1 2001/03/15 19:02:08 radphi Exp $";

int main(int nargs, char *args[]){
  int length;
  int run_number;
  int i;
  unsigned char *char_array=NULL;
  if(nargs!=6){
    fprintf(stderr,"usage:\n");
    fprintf(stderr,"map_fill_char_item_with_zeros <map> <subsys> <item> <length> <runNo>\n");
    exit(1);
  }
  sscanf(args[4],"%d",&length);
  sscanf(args[5],"%d",&run_number);
  char_array=malloc(sizeof(unsigned char)*length);
  if(char_array==NULL){
    fprintf(stderr,"There was a problem with the malloc.\n");
    exit(1);
  }
  for(i=0;i<length;i++){
    char_array[i]=0;
  }
  if(map_put_char(args[1],args[2],args[3],length,char_array,run_number)){
    fprintf(stderr,"map_fill_char_item_with_zeros:\n");
    fprintf(stderr,"encountered error writing into map.\n");
    exit(1);
  }
  fprintf(stderr,"Successfully wrote initialization into map!\n");
  
}
