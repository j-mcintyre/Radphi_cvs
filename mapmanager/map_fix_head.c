/* 
   map_fix_head.c
   Created           : 19-NOV-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_FIX_HEAD                                                      */
/*   -----------                                                       */
/*                                                                     */
/*         Created     : 19-NOV-1993    Author : Thom Sulanke          */
/*         Purpose     : fix the subsystem table header from version 1 */
/*                       to version 2                                  */
/*                                                                     */
/***********************************************************************/
 
int map_fix_head( const char filename[], int nsubsystems)
     
{
  int fd ;
  pointer_t subsystem_loc, header_loc;
  int status;
  nametable_t *table;
  subsystem_t subsystem;
  subtableheader_t header;
  
  /* open file for read and write */
  
  fd = map_open_rw(filename);
  if ( fd < 0 )
    {
      map_syserror(fd,
		   "can't open TheMap (%s) for fixing subsystem",
		   filename);
      return fd;
    }
  
  /* find location of first subsystem */

  status = map_read_pointer(&subsystem_loc,sizeof header_loc,fd,ROOT);
  if ( status != MAP_OK )
    {
      map_close_map(fd);
      return status;
    }
  
  if ( subsystem_loc != NULL_LOC )
    {
      status = map_read_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
    }
  else
    subsystem.next = NULL_LOC;
  
  if ( subsystem.next < 0 )
    /* really the subsystem table header */
    printf("Subsystem table header already fixed.\n");
  else
    {
      if ( nsubsystems < INIT_NAMETABLE_SIZE )
	header.table_length = INIT_NAMETABLE_SIZE;
      else
	header.table_length = nsubsystems;
      table = (nametable_t *) calloc(header.table_length,
				     sizeof(nametable_t));
      header.table_used = 0;
      while ( subsystem_loc != NULL_LOC )
	{
	  status = map_read_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	  memcpy(table[header.table_used].name,subsystem.name,MAP_NAMELENGTH);
	  table[header.table_used].loc = subsystem_loc;
	  header.table_used++;
	  subsystem_loc = subsystem.next;
	}  /* While */
      
      status = map_write_nametable(table,sizeof(nametable_t)*header.table_length,
			           fd,&header.table);
      free(table);
      if ( status < 0 )
	{
	  map_close_map(fd);
	  return status;
	}
      header.negversion = -VERSION;
      header.nullname = 0;
      status = map_write_subtableheader(&header,sizeof header,fd,&header_loc);
      if ( status < MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
      status = map_overwrite_pointer(&header_loc,sizeof header_loc,fd,ROOT);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
    }
           
  return map_close_map(fd);
  
}
