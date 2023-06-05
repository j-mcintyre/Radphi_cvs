/* 
   map_rem_item.c
   Created           :  30-Sep-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_REM_ITEM                                                      */
/*   ------------                                                      */
/*                                                                     */
/*         Created     :  30-SEP-1993    Author : Thom Sulanke         */
/*         Purpose     :  remove an item from map                      */
/*                                                                     */
/***********************************************************************/
 
int map_rem_item( const char filename[], const char subsystemname[], 
		 const char itemname[] )
{
  int fd ;
  pointer_t pred_loc, item_loc, subsystem_loc, succ_loc, header_loc;
  item_t item;
  int status, tindex;
  subtableheader_t header;
  subsystem_t subsystem;
  nametable_t *nametable;
  int ret;
  
  /* check arguments */
  
  if ( strlen(subsystemname) > MAP_NAMELENGTH-1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_rem_item: Subsystem name (%s) longer than %d characters.",
		subsystemname, MAP_NAMELENGTH-1);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  if ( strlen(itemname) > MAP_NAMELENGTH-1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_rem_item: Item name (%s) longer than %d characters.",
		itemname, MAP_NAMELENGTH-1);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  /* open file for read and write */
  
  fd = map_open_rw(filename);
  if ( fd < MAP_OK )
    {
      map_syserror(fd,
		   "map_rem_item: can't open (for read/write) TheMap (%s).",
		   filename);
      return fd;
    }
  
  /* find location of item in list */
  
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
    ret = map_find_same(fd,ROOT,subsystemname,0,NULL_LOC,&tindex,&succ_loc,&subsystem_loc);
  else
    ret = map_find_same(fd,NULL_LOC,subsystemname,
				  header.table_used,header.table,&tindex,&succ_loc,
				  &subsystem_loc);
  
  if ( ret < 0 )
    {
      if ( ret == MAP_USER_WARN_NOMATCH )
	map_warn(ret,
		 "map_rem_item: subsystem (%s) does not exist",
		 subsystemname);
      map_close_map(fd);
      return ret;
    }
  
  status = map_read_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
  if ( status != MAP_OK )
    {
      map_close_map(fd);
      return status;
    }
  
  if ( subsystem.list >= 0 )

    /* version 1 */

    {
      ret = map_find_pred(fd,subsystem_loc+SUBSYSTEM_LIST_OFFSET,
			       itemname,0,NULL_LOC,&tindex,&succ_loc,&pred_loc);
      if ( ret < 0 )
	{
	  map_close_map(fd);
	  return ret;
	}
      status = map_read_pointer(&item_loc,sizeof item_loc,fd,pred_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
      if ( item_loc != NULL_LOC )
	{
	  status = map_read_item(&item,sizeof item,fd,item_loc);
	  if ( status != MAP_OK )
	    {
	      map_close_map(fd);
	      return status;
	    }
	  if ( strncmp(itemname, item.name, MAP_NAMELENGTH) != 0 )
	    {
	      map_warn(MAP_USER_WARN_NOMATCH,
		       "map_rem_item: item (%s) does not exists for subsystem (%s)",
		       itemname,
		       subsystemname);
	      map_close_map(fd);
	      return MAP_USER_WARN_NOMATCH;
	    }
	}
      else
	{
	  map_warn(MAP_USER_WARN_NOMATCH,
		   "map_rem_item: item (%s) does not exists for subsystem (%s)",
		   itemname,
		   subsystemname);
	  map_close_map(fd);
	  return MAP_USER_WARN_NOMATCH;
	}
      
      /* unlink entry */
      
      /* update pointer to new entry */
      
      status = map_overwrite_pointer(&item.next,sizeof item.next,fd,pred_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}

  
    }
  else
    
    /* version 2 */
    
    {
      ret = map_find_same(fd,NULL_LOC,itemname,
			       subsystem.table_used,subsystem.table,&tindex,
			       &succ_loc,&item_loc);
      if ( ret < 0 )
	{
	  if ( ret == MAP_USER_WARN_NOMATCH )
	    map_warn(ret,
		 "map_rem_item: item (%s) does not exist in subsystem (%s)",
		     itemname,
		     subsystemname);
	  map_close_map(fd);
	  return ret;
	}
      if ( subsystem.table_used-tindex-1 > 0 )
	{
	  /* read portion of table after entry and rewrite it over entry */

	  nametable = (nametable_t *) calloc(subsystem.table_used-tindex-1,
					     sizeof(nametable_t));
	  status = map_read_nametable(nametable,
			              (subsystem.table_used-tindex-1)*sizeof(nametable_t),
			              fd,subsystem.table+(tindex+1)*(int)sizeof(nametable_t));
	  if ( status != MAP_OK )
	    {
	      free(nametable);
	      map_close_map(fd);
	      return status;
	    }
	  status = map_overwrite_nametable(nametable,
				 (subsystem.table_used-tindex-1)*sizeof(nametable_t),
				 fd,
				 subsystem.table+tindex*(int)sizeof(nametable_t));
	  if ( status != MAP_OK )
	    {
	      free(nametable);
	      map_close_map(fd);
	      return status;
	    }
	  free(nametable);
	}
      subsystem.table_used--;
      status = map_overwrite_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}

    } 

  /* close the file */
  
  status = map_close_map(fd);
  if ( status != MAP_OK ) return status;
  
  /* if the item removed was the last of this subsystem, remove the subsystem as well */

  if ( subsystem.list >= 0 )
    {
      if ( (pred_loc == subsystem_loc+SUBSYSTEM_LIST_OFFSET) &&
	  (item.next == NULL_LOC) )
	return map_rem_sub(filename,subsystemname);
    }
  else
    {
      if ( subsystem.table_used == 0 )
	return map_rem_sub(filename,subsystemname);

    }

  return MAP_OK;
    
}

