/* 
   mapfixit.c
   Created           : 16-NOV-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <limits.h>
#include <map_manager.h>
#include <map_internal.h>
#include <stdlib.h>

static char  rcsid[] = "$Id: mapfixit.c,v 1.1 1996/11/08 15:14:05 lfcrob Exp $";

/* $Log: mapfixit.c,v $
/* Revision 1.1  1996/11/08 15:14:05  lfcrob
/* Initial revision
/* */

int	map_fix_item(const char filename[], const char subsystemname[], 
	    const char itemname[], int narray);
char *progname;
 

/***********************************************************************/
/*                                                                     */
/*   mapfixit                                                          */
/*   -------                                                           */
/*                                                                     */
/*         Created     : 16-NOV-1993    Author : Thom Sulanke          */
/*         Purpose     : convert map from version 1 to version 2       */
/*                                                                     */
/*                                                                     */
/***********************************************************************/
 
main( argc, argv )
     int argc;
     char *argv[];
{
  char subsystemname[MAP_NAMELENGTH], itemname[MAP_NAMELENGTH];
  char nextsubsystemname[MAP_NAMELENGTH], nextitemname[MAP_NAMELENGTH];
  item_t item;
  int i;                /* index */
  int atime, vtime;
  int length;
  int type;
  int nvals;
  int nitems;
  int nsubsystems;
  
  progname = argv[0];
  
  if ( argc != 2 )
    {
      printf("Usage: %s themap\n", progname);
      exit(0);
    }
  
  printf("Fixing: %s.\n",argv[1]);

  /* add subsystem table header */

  map_stat_head(argv[1], &nsubsystems);
  printf("Number of subsytems: %d.\n",nsubsystems);

  /* loop through all subsystems and items */
  
  for ( map_next_sub( argv[1], "*", nextsubsystemname, &nitems) ; 
       strcmp(nextsubsystemname,"*") !=0 ; 
       map_next_sub(argv[1], subsystemname, nextsubsystemname, &nitems) )
    {
      memcpy(subsystemname, nextsubsystemname, MAP_NAMELENGTH);
      printf("Subsystem: %s.  Number of items: %d.\n",subsystemname,nitems);
      for ( map_next_item( argv[1], subsystemname, "*", nextitemname, 
			  &length, &type, &nvals) ;
	   strcmp(nextitemname,"*") !=0 ; 
	   map_next_item(argv[1], subsystemname, itemname, nextitemname, 
			 &length, &type, &nvals))
        {
	  memcpy(itemname, nextitemname, MAP_NAMELENGTH);
	  printf("  Item: %s.  Array size: %d.  Array type: %d.  Number of arrays: %d.\n",itemname,length,type,nvals);
	  map_fix_item(argv[1], subsystemname, itemname, nvals);

	}  /* For */
      
      map_fix_sub(argv[1], subsystemname, nitems);
    
    }  /* For */
  
  map_fix_head(argv[1], nsubsystems);
  
  return 0;
}
