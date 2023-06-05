/* 
   map_stat_head.c
   Created           : 19-NOV-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include "map_manager.h"
#include "map_internal.h"
 

/***********************************************************************/
/*                                                                     */
/*   MAP_STAT_HEAD                                                     */
/*   -------------                                                     */
/*                                                                     */
/*         Created     : 19-NOV-1993    Author : Thom Sulanke          */
/*         Purpose     : return number of subsystems                   */
/*                                                                     */
/***********************************************************************/
 
int map_stat_head( const char filename[], int *nsubsystems)
     
{
  int fd ;
  pointer_t pred_loc, header_loc;
  int status;
  subtableheader_t header;
    
  /* set default for early return */
  
  *nsubsystems = 0;
  
  /* open file for read */
  
  fd = map_open_ro(filename);
  if ( fd < 0 )
    {
      map_syserror(fd,
		   "can't open TheMap (%s) for finding next subsystem",
		   filename);
      return fd;
    }
    
  status = map_read_pointer(&header_loc,sizeof header_loc,fd,ROOT);
  if ( status != MAP_OK )
    {
      map_close_map(fd);
      return status;
    }
  
  status = map_read_subtableheader(&header,sizeof header,fd,header_loc);
  if ( status != MAP_OK )
    {
      map_close_map(fd);
      return status;
    }
  
  if ( header.negversion >= 0 )
    
    /* version 1 */
    
    {
      *nsubsystems = 0;
      pred_loc = header_loc;
      while ( pred_loc != NULL_LOC )
	{
	  (*nsubsystems)++;
	  status = map_read_pointer(&pred_loc,sizeof pred_loc,fd,pred_loc);
	  if ( status != MAP_OK ) 
	    {
	      map_close_map(fd);
	      return status;
	    }
	}
    }
  else

    /* version 2 */
    
    *nsubsystems = header.table_used;

  return map_close_map(fd);
  
}
