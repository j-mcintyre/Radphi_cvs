/* 
   map_rem_arr.c
   Created           : 29-NOV-1993 by Thom Sulanke
*/

#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>

/***********************************************************************/
/*                                                                     */
/*   MAP_REM_ARR                                                       */
/*   -----------                                                       */
/*                                                                     */
/*         Created     : 29-NOV-1993    Author : Thom Sulanke          */
/*         Purpose     : remove an array of values from database.      */
/*                                                                     */
/*         Input:    filename: the map file name                       */
/*                     subsystemname: subsystem name                   */
/*                     itemname: item name                             */
/*                     atime: first time when values to be removed     */
/*                            were valid                               */
/***********************************************************************/
 
int map_rem_arr( const char filename[], const char subsystemname[], 
		const char itemname[], int atime)
     
{
  int fd ;
  pointer_t item_loc, pred_loc, suc_loc;
  item_t item;
  int status, tindex, suc_time;
  table_t *table;
  
  /* check arguements */
  
  if ( strlen(subsystemname) > MAP_NAMELENGTH-1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_rem_arr: Subsystem name (%s) longer than %d characters.",
		subsystemname,
		MAP_NAMELENGTH-1);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  if ( strlen(itemname) > MAP_NAMELENGTH-1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_rem_arr: Item name (%s) longer than %d characters.",
		itemname,
		MAP_NAMELENGTH-1);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  /* open file for read and write */
  
  fd = map_open_rw(filename);
  if ( fd < MAP_OK )
    {
      map_syserror(fd,
		   "map_rem_arr: can't open TheMap (%s) for adding values",
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
  
  /* find array in list */
  
  status = map_time_pred(fd,item_loc+ITEM_LIST_OFFSET,item,atime,
			   &tindex,&suc_loc,&suc_time,&pred_loc);
  if ( status < 0 )
    {
      map_close_map(fd);
      return status;
    }
  
  /* check if sucessor is same time */
  
  if ( suc_loc == NULL_LOC || suc_time != atime )
    {
      map_warn(MAP_USER_WARN_NOFIND,
		"map_rem_arr: No values exist for this time (%d) \n for subsystem (%s), item (%s).",
		atime,
		subsystemname,
		itemname);
	  map_close_map(fd);
	  return MAP_USER_WARN_NOFIND;
    }
      
  /* update table */
    
  if ( tindex+1 != item.table_used-1 )

    /* if array is not last on list, read end of table and rewrite it */

    {
      table = (table_t *) calloc(item.table_used-tindex-2,sizeof(table_t));
      status = map_read_table(table,(item.table_used-tindex-2)*sizeof(table_t),
		 	      fd,item.table+(tindex+2)*(int)sizeof(table_t));
      if ( status != MAP_OK )
	{
	  free(table);
	  map_close_map(fd);
	  return status;
	}
      status = map_overwrite_table(table,
			     (item.table_used-tindex-2)*sizeof(table_t),
			     fd,item.table+(tindex+1)*(int)sizeof(table_t));
      if ( status != MAP_OK )
	{
	  free(table);
	  map_close_map(fd);
	  return status;
	}
      free(table);
    }
  item.table_used--;
  status = map_overwrite_item(&item,sizeof item,fd,item_loc);
  if ( status != MAP_OK )
    {
      map_close_map(fd);
      return status;
    }
  
  /* close the file */
  
  return map_close_map(fd);
  
}
