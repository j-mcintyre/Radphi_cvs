/* 
   map_put_float.c
   Created           : 12-NOV-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>

/***********************************************************************/
/*                                                                     */
/*   MAP_PUT_FLOAT                                                     */
/*   -------------                                                     */
/*                                                                     */
/*         Created     : 12-NOV-1992    Author : Thom Sulanke          */
/*         Purpose     : put array of float values into the map.       */
/*                                                                     */
/***********************************************************************/
 
int map_put_float( const char filename[], const char subsystemname[], 
		  const char itemname[], int length, const float array[], 
		  int atime )
 
{
  int fd ;
  pointer_t item_loc, pred_loc, array_loc;
  item_t item;
  arrayheader_t header;
  int status, tindex, suc_time;
  table_t *table;
  pointer_t dummy;
  
  header.spare[0] = 0;
  
  /* get the current time */
  
  if ( atime == 0 )
    {
      header.time = time(NULL);
    }
  else
    {
      header.time = abs(atime);
    }
  
  /* check arguements */
  
  if ( strlen(subsystemname) > MAP_NAMELENGTH-1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_put_float: Subsystem name (%s) longer than %d characters.",
		subsystemname,
		MAP_NAMELENGTH-1);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  if ( strlen(itemname) > MAP_NAMELENGTH-1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_put_float: Item name (%s) longer than %d characters.",
		itemname,
		MAP_NAMELENGTH-1);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  /* open file for read and write */
  
  fd = map_open_rw(filename);
  if ( fd < MAP_OK )
    {
      map_syserror(fd,
		   "map_put_float: can't open TheMap (%s) for adding values",
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
  
  /* check array size */
  
  if ( length != item.length )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_put_float: Array length (%d) incorrect for subsystem (%s), item (%s).\n Should be %d.",
		length,
		subsystemname,
		itemname,
		item.length);
      map_close_map(fd);
      return MAP_USER_ERROR_ARGUMENT;
    }
  header.length = length;
  
  /* check array type */
  
  if ( item.type != 1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_put_float: Array type (%d) for subsystem (%s), item (%s) \n requires call to different map_put_* routine.",
		item.type,
		subsystemname,
		itemname);
      map_close_map(fd);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  /* find place for array in list */
  
  status = map_time_pred(fd,item_loc+ITEM_LIST_OFFSET,item,header.time,
			   &tindex,&header.next,&suc_time,&pred_loc);
  if ( status < 0 )
    {
      map_close_map(fd);
      return status;
    }
  
  /* check if sucessor is same time */
  
  if ( suc_time == header.time )
    {
      if ( atime >= 0 )
	{
	  map_warn(MAP_USER_WARN_NOREPLACE,
		    "map_put_float: Values already exist for this time (%d) \n for subsystem (%s), item (%s).",
		    header.time,
		    subsystemname,
		    itemname);
	  map_close_map(fd);
	  return MAP_USER_WARN_NOREPLACE;
	}
      else
	{
	  
	  status = map_overwrite_float(array,item.length*sizeof(float),fd,header.next+ARRAY_VALUES_OFFSET);
	  if ( status != MAP_OK )
	    {
	      map_close_map(fd);
	      return status;
	    }
	}
      
    }
  else
    {
      
      /* fill in entry next pointer */
      
      status = map_write_arrayheader(&header,sizeof header,fd,&array_loc);
      if ( status < 0 )
	{
	  map_close_map(fd);
	  return status;
	}
      
      /* insert entry */
      
      status = map_write_float(array,item.length*sizeof(float),fd,&dummy);
      if ( status < 0 )
	{
	  map_close_map(fd);
	  return status;
	}
      
      if ( item.table_length > 0 )
	
	/* version 2 format */
	/* update table */
	  
	{
	  if ( item.table_used < item.table_length )
	    table = (table_t *) calloc(item.table_used+1,sizeof(table_t));
	  else
	    table = (table_t *) calloc(item.table_length * TABLE_FACTOR,sizeof(table_t));
	  
	  if ( item.table_used > 0 )
	    {
	      status = map_read_table(table,item.table_used*sizeof(table_t),fd,item.table);
	      if ( status != MAP_OK )
		{
		  free(table);
		  map_close_map(fd);
		  return status;
		}
	    }
	  if ( item.table_used-tindex-1 > 0 )
	    memmove(&table[tindex+2],&table[tindex+1],(item.table_used-tindex-1)*2*sizeof(int));
	  table[tindex+1].time = header.time;
	  table[tindex+1].loc = array_loc;
	  item.table_used++;
	  if ( item.table_used <= item.table_length )
	    {
	      status = map_overwrite_table(table,item.table_used*sizeof(table_t),fd,item.table);
	      if ( status != MAP_OK )
		{
		  free(table);
		  map_close_map(fd);
		  return status;
		}
	    }
	  else
	    {
	      item.table_length = item.table_length * TABLE_FACTOR;
	      status = map_write_table(table,item.table_length*sizeof(table_t),
				       fd,&item.table);
	      if ( status < 0 )
		{
		  free(table);
		  map_close_map(fd);
		  return status;
		}
	    }
	  free(table);
	  status = map_overwrite_item(&item,sizeof item,fd,item_loc);
	  if ( status != MAP_OK )
	    {
	      map_close_map(fd);
	      return status;
	    }
	}
 
      else
	
	/* version 1 format */
	/* update pointer to new entry */
	
	{
	  status = map_overwrite_pointer(&array_loc,sizeof array_loc,fd,pred_loc);
	  if ( status != MAP_OK )
	    {
	      map_close_map(fd);
	      return status;
	    }
	  
	}
      
    }
  
  /* close the file */
  
  return map_close_map(fd);
  
}
