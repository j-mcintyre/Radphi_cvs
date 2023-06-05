/* 
   map_dump.c
   Created           : 18-DEC-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   map_dump                                                          */
/*   --------                                                          */
/*                                                                     */
/*         Created     : 18-DEC-1992    Author : Thom Sulanke          */
/*         Purpose     : dump the map                                  */
/*                                                                     */
/***********************************************************************/
 
int map_dump( const char filename[] )
{
  char subsystemname[MAP_NAMELENGTH], itemname[MAP_NAMELENGTH];
  char nextsubsystemname[MAP_NAMELENGTH], nextitemname[MAP_NAMELENGTH];
  int i;                /* index */
  int *array;
  char *carray;
  float *farray;
  int atime, vtime;
  time_t vtime_time;
  int length,olength;
  int type;
  int nvals;
  int nitems;
  int nsubsystems;

  /* get number of subsystems */

  map_stat_head( filename, &nsubsystems );
  printf("Map: %s.  Number of subsystems: %d.\n",filename,nsubsystems);

  /* dump list of subsystems and items */
  
  for ( map_next_sub( filename, "*", nextsubsystemname, &nitems) ; 
       strcmp(nextsubsystemname,"*") !=0 ; 
       map_next_sub(filename, subsystemname, nextsubsystemname, &nitems) )
    {
      memcpy(subsystemname, nextsubsystemname, MAP_NAMELENGTH);
      printf("Subsystem: %s.  Number of items: %d.\n",subsystemname,nitems);
      for ( map_next_item( filename, subsystemname, "*", nextitemname, 
			  &length, &type, &nvals) ;
	   strcmp(nextitemname,"*") !=0 ; 
	   map_next_item(filename, subsystemname, itemname, nextitemname, 
			 &length, &type, &nvals))
        {
	  
	  char typeStr[12];
	  memcpy(itemname, nextitemname, MAP_NAMELENGTH);
	  switch(type){
	  case MAP_INT:
	    sprintf(typeStr,"int");
	    break;
	  case MAP_FLOAT:
	    sprintf(typeStr,"float");
	    break;
	  case MAP_CHAR:
	    sprintf(typeStr,"char");
	    break;
	  default:
	    sprintf(typeStr,"unknown");
	    break;
	  }
	  
	  printf("  Item: %s.  Array size: %d.  Array type: %s.  Number of arrays: %d.\n",itemname,length,typeStr,nvals);
	  atime = INT_MAX;
	  switch ( type ) {
	  case 0: 
	    array = (int *) calloc(length,sizeof(int));
	    if ( array == NULL )
	      {
		printf("map_dump: error allocating memory for array.\n");
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_int( filename, subsystemname, itemname, length,
			    array, atime, &vtime);
		vtime_time = vtime;
		printf("    Time: %9d, %.24s, first value: %d, last value: %d.\n",
		       vtime, ctime(&vtime_time), array[0], array[length-1]);
		atime = vtime - 1;
	      }  /* For */
	    free(array);      
	    break;
	  case 1: 
	    farray = (float *) calloc(length,sizeof(float));
	    if ( farray == NULL )
	      {
		printf("map_dump: error allocating memory for array.\n");
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_float( filename, subsystemname, itemname, length,
			      farray, atime, &vtime);
		vtime_time = vtime;
		printf("    Time: %9d, %.24s, first value: %f, last value: %f.\n",
		       vtime, ctime(&vtime_time), farray[0], farray[length-1]);
		atime = vtime - 1;
	      }  /* For */
	    free(farray);      
	    break;
	  case 2: 
	    carray = (char *) calloc(length,sizeof(char));
	    if ( carray == NULL )
	      {
		printf("map_dump: error allocating memory for array.\n");
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_char( filename, subsystemname, itemname, length,
			     carray, atime, &vtime, &olength);
		vtime_time = vtime;
		printf("    Time: %9d, %.24s, Length: %5d, first values in Hex: ",
		       vtime,ctime(&vtime_time),olength);
		if (olength > 10) olength = 10;
		for(i=0;i<olength;i++) printf(" %x",carray[i]);
		printf("\n");
		atime = vtime - 1;
	      }  /* While */
	    free(carray);      
	    break;
	  }
        }  /* For */
    }  /* For */
  
  return 0;
}
