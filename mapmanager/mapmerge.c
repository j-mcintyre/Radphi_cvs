/* 
        mapmerge.c
         Created           : 15-NOV-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <limits.h>
#include <map_manager.h>
#include <map_internal.h>
#include <stdlib.h>

static char  rcsid[] = "$Id: mapmerge.c,v 1.1 1996/11/08 15:14:09 lfcrob Exp $";

/* $Log: mapmerge.c,v $
/* Revision 1.1  1996/11/08 15:14:09  lfcrob
/* Initial revision
/* */

char *progname;
 

/***********************************************************************/
/*                                                                     */
/*   mapmerge                                                          */
/*   --------                                                          */
/*                                                                     */
/*         Created     : 15-NOV-1993    Author : Thom Sulanke          */
/*         Purpose     : merge one map into another                    */
/*                                                                     */
/***********************************************************************/
 
main( argc, argv )
     int argc;
     char *argv[];
{
    char subsystemname[MAP_NAMELENGTH], itemname[MAP_NAMELENGTH];
    char nextsubsystemname[MAP_NAMELENGTH], nextitemname[MAP_NAMELENGTH];
    int i;                /* index */
    int *array;
    char *carray;
    float *farray;
    int atime, vtime;
    int length,olength;
    int type;
    int nvals;
    int nitems;

    progname = argv[0];

    if ( argc != 3 )
    {
        printf("Usage: %s frommap intomap\n", progname);
	exit(0);
    }

/* loop through each subsystem of frommap  */
 
    for ( map_next_sub( argv[1], "*", nextsubsystemname, &nitems) ; 
        strcmp(nextsubsystemname,"*") !=0 ; 
        map_next_sub(argv[1], subsystemname, nextsubsystemname, &nitems) )
    {
      memcpy(subsystemname, nextsubsystemname, MAP_NAMELENGTH);

      fprintf(stderr,"Subsystem: %s\r",subsystemname);

      for ( map_next_item( argv[1], subsystemname, "*", nextitemname, 
			  &length, &type, &nvals) ;
	   strcmp(nextitemname,"*") !=0 ; 
	   map_next_item(argv[1], subsystemname, itemname, nextitemname, &length, &type, &nvals))
        {
	  memcpy(itemname, nextitemname, MAP_NAMELENGTH);

	  fprintf(stderr,"Subsystem: %s, item: %s                                \r",subsystemname,itemname);

	  map_add_item1(argv[2], subsystemname, itemname, length, type, nvals);

	  atime = INT_MAX;
	  switch ( type ) {
	  case 0: 
	    array = (int *) calloc(length,sizeof(int));
	    if ( array == NULL )
	      {
		printf("%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_int( argv[1], subsystemname, itemname, length,
			    array, atime, &vtime);
		map_put_int( argv[2], subsystemname, itemname, length,
			    array, -vtime);
		atime = vtime - 1;
	      }  /* For */
	    free(array);      
	    break;
	  case 1: 
	    farray = (float *) calloc(length,sizeof(float));
	    if ( farray == NULL )
	      {
		printf("%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_float( argv[1], subsystemname, itemname, length,
			      farray, atime, &vtime);
		map_put_float( argv[2], subsystemname, itemname, length,
			      farray, -vtime);
		atime = vtime - 1;
	      }  /* For */
	    free(farray);      
	    break;
	  case 2: 
	    carray = (char *) calloc(length,sizeof(char));
	    if ( carray == NULL )
	      {
		printf("%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_char( argv[1], subsystemname, itemname, length,
			     carray, atime, &vtime, &olength);
/*
printf("nvals: %d, vtime: %d, olength: %d\n",nvals,vtime,olength);
*/
		map_put_char( argv[2], subsystemname, itemname, olength,
			     carray, -vtime);
		atime = vtime - 1;
	      }  /* While */
	    free(carray);      
	    break;
	  }
        }  /* For */
    }  /* For */
 
    return 0;
}
