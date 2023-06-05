/* 
   map_fix_sub.c
   Created           : 19-NOV-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_FIX_SUB                                                       */
/*   -----------                                                       */
/*                                                                     */
/*         Created     : 19-NOV-1993    Author : Thom Sulanke          */
/*         Purpose     : fix one subsystem from version 1 to version 2 */
/*                                                                     */
/***********************************************************************/
 
int map_fix_sub( const char filename[], const char subsystemname[], 
		 int nitems)
     
{
  int fd ;
  pointer_t subsystem_loc, item_loc, table_loc;
  int status;
  nametable_t *table;
  item_t item;
  subsystem_t subsystem;
  
  /* open file for read and write */
  
  fd = map_open_rw(filename);
  if ( fd < 0 )
    {
      map_syserror(fd,
		   "can't open TheMap (%s) for fixing subsystem",
		   filename);
      return fd;
    }
  
  /* find subsystem location and read subsystem */
  
  status = map_find_sub(fd,subsystemname,&subsystem_loc);
  if ( status < 0 )
    {
      map_close_map(fd);
      return status;
    }
  status = map_read_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
  if ( status != MAP_OK ) 
    {
      map_close_map(fd);
      return status;
    }
  
  /* check if already fixed */
  
  if ( subsystem.list == -VERSION )
    printf("    Subsystem already fixed.\n");
  else
    {
      if ( nitems < INIT_NAMETABLE_SIZE )
	subsystem.table_length = INIT_NAMETABLE_SIZE;
      else
	subsystem.table_length = nitems;
      table = (nametable_t *) calloc(subsystem.table_length,
				     sizeof(nametable_t));
      subsystem.table_used = 0;
      item_loc = subsystem.list;
      while ( item_loc != NULL_LOC )
	{
	  status = map_read_item(&item,sizeof item,fd,item_loc);
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	  memcpy(table[subsystem.table_used].name,item.name,MAP_NAMELENGTH);
	  table[subsystem.table_used].loc = item_loc;
	  subsystem.table_used++;
	  item_loc = item.next;
	}  /* While */
      
      status = map_write_nametable(table,sizeof(nametable_t)*subsystem.table_length,
			           fd,&table_loc);
      free(table);
      if ( status < 0 )
	{
	  map_close_map(fd);
	  return status;
	}
      subsystem.table = table_loc;
      subsystem.list = -VERSION;
      status = map_overwrite_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
    }
           
  return map_close_map(fd);
  
}
