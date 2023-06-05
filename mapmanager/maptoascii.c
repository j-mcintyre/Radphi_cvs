/* 
        maptoascii.c
         Created           : 11-May-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <map_manager.h>
#include <map_internal.h>
#include <stdlib.h>

static char  rcsid[] = "$Id: maptoascii.c,v 1.2 2005/03/30 12:40:48 radphi Exp $";

/* $Log: maptoascii.c,v $
/* Revision 1.2  2005/03/30 12:40:48  radphi
/* maptoascii.c
/*    - modified format specifier %f to %g to prevent loss of precision
/*      in exporting a map file in case of constants of order < 1e5  [rtj]
/*
/* Revision 1.1  1996/11/08 15:14:26  lfcrob
/* Initial revision
/* */


char *progname;
 

/***********************************************************************/
/*                                                                     */
/*   maptoascii                                                        */
/*   ----------                                                        */
/*                                                                     */
/*         Created     : 11-May-1993    Author : Thom Sulanke          */
/*         Purpose     : write information in ascii for transfer to    */
/*                       incompatible systems                          */
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
    int length, olength;
    int type;
    int nvals;
    int nitems;

    progname = argv[0];

    if ( argc != 2 )
    {
        fprintf(stderr,"Usage: %s themap\n", progname);
	exit(0);
    }

/* write list of subsystems and items */
 
    for ( map_next_sub( argv[1], "*", nextsubsystemname, &nitems) ; 
        strcmp(nextsubsystemname,"*") !=0 ; 
        map_next_sub(argv[1], subsystemname, nextsubsystemname, &nitems) )
    {
      memcpy(subsystemname, nextsubsystemname, MAP_NAMELENGTH);
      printf("S ");
      for(i=0;i<MAP_NAMELENGTH;i++){printf(" %x",subsystemname[i]);}
	  printf("\n");
      for ( map_next_item( argv[1], subsystemname, "*", nextitemname, &length, &type, &nvals) ;
	   strcmp(nextitemname,"*") !=0 ; 
	   map_next_item(argv[1], subsystemname, itemname, nextitemname, &length, &type, &nvals))
        {
	  memcpy(itemname, nextitemname, MAP_NAMELENGTH);
	  printf("I %d %d",length,type);
      for(i=0;i<MAP_NAMELENGTH;i++){printf(" %x",itemname[i]);}
	  printf("\n");
	  atime = 0;
	  switch ( type ) {
	  case 0: 
	    array = (int *) calloc(length,sizeof(int));
	    if ( array == NULL )
	      {
		fprintf(stderr,"%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_int( argv[1], subsystemname, itemname, length,
			    array, atime, &vtime);
		printf("A %d", vtime);
		for(i=0;i<length;i++){
		  if ( i % 10 == 0 )
		    printf("\n");
		  printf(" %d",array[i]);
		}
		printf("\n");
		atime = vtime - 1;
	      }  /* For */
	    free(array);      
	    break;
	  case 1: 
	    farray = (float *) calloc(length,sizeof(float));
	    if ( farray == NULL )
	      {
		fprintf(stderr,"%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_float( argv[1], subsystemname, itemname, length,
			      farray, atime, &vtime);
		printf("F %d", vtime);
		for(i=0;i<length;i++){
		  if ( i % 10 == 0 )
		    printf("\n");
		  printf(" %g",farray[i]);
		}
		printf("\n");
		atime = vtime - 1;
	      }  /* For */
	    free(farray);      
	    break;
	  case 2: 
	    carray = (char *) calloc(length,sizeof(char));
	    if ( carray == NULL )
	      {
		fprintf(stderr,"%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    for ( ; nvals>0 ; nvals-- )
	      {
		map_get_char( argv[1], subsystemname, itemname, length,
			     carray, atime, &vtime, &olength);
		printf("C %d %d", vtime, olength);
		for(i=0;i<olength;i++){
		  if ( i % 20 == 0 )
		    printf("\n");
		  printf(" %x",carray[i]);
		}
		printf("\n");
		atime = vtime - 1;
	      }  /* For */
	    free(carray);      
	    break;
	  }
	  printf("E\n"); /* end of values */
        }  /* For */
      printf("E\n"); /* end of items */
    }  /* For */
    printf("E\n"); /* end of subsystems */
    
    return 0;
}
