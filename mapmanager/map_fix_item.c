/* 
   map_fix_item.c
   Created           : 16-NOV-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <string.h>
#include <map_manager.h>
#include <map_internal.h>
#include <stdlib.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_FIX_ITEM                                                      */
/*   ------------                                                      */
/*                                                                     */
/*         Created     : 16-NOV-1993    Author : Thom Sulanke          */
/*         Purpose     : fix one item from version 1 to version 2      */
/*                                                                     */
/***********************************************************************/
 
int map_fix_item( const char filename[], const char subsystemname[], 
		 const char itemname[], int nvals)
     
{
  int fd ;
  pointer_t item_loc, entry_loc, table_loc;
  int status;
  item_t item;
  table_t *table;
  arrayheader_t entry;
  
  /* open file for read and write */
  
  fd = map_open_rw(filename);
  if ( fd < 0 )
    {
      map_syserror(fd,
		   "can't open TheMap (%s) for fixing item",
		   filename);
      return fd;
    }
  
  /* find item location and read item */
  
  status = map_find_item(fd,subsystemname,itemname,&item_loc);
  if ( status < 0 )
    {
      map_close_map(fd);
      return status;
    }
  status = map_read_item(&item,sizeof item,fd,item_loc);
  if ( status != MAP_OK ) 
    {
      map_close_map(fd);
      return status;
    }
  
  /* check if already fixed */
  
  if ( item.list == -VERSION )
    printf("    Item already fixed.\n");
  else
    {
      if ( nvals < INIT_TABLE_SIZE )
	item.table_length = INIT_TABLE_SIZE;
      else
	item.table_length = nvals;
      table = (table_t *) calloc(item.table_length,sizeof(table_t));
      item.table_used = 0;
      entry_loc = item.list;
      while ( entry_loc != NULL_LOC )
	{
	  status = map_read_arrayheader(&entry,sizeof entry,fd,entry_loc);
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	  table[item.table_used].time = entry.time;
	  table[item.table_used].loc = entry_loc;
	  item.table_used++;
	  entry_loc = entry.next;
	}  /* While */
      
      status = map_write_table(table,sizeof(table_t)*item.table_length,fd,&table_loc);
      free(table);
      if ( status < 0 )
	{
	  map_close_map(fd);
	  return status;
	}
      item.table = table_loc;
      item.list = -VERSION;
      status = map_overwrite_item(&item,sizeof item,fd,item_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
    }
           
  return map_close_map(fd);
  
}
