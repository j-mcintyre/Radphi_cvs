/* 
   map_rem_sub.c
   Created           :  30-Sep-1993 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_REM_SUB                                                       */
/*   ------------                                                      */
/*                                                                     */
/*         Created     :  30-SEP-1993    Author : Thom Sulanke         */
/*         Purpose     :  remove a subsystem from map                  */
/*                                                                     */
/***********************************************************************/
 
int map_rem_sub( const char filename[], const char subsystemname[] )
{
  int fd ;
  pointer_t pred_loc, subsystem_loc, succ_loc, header_loc;
  subsystem_t subsystem;
  int status, tindex;
  nametable_t *table;
  subtableheader_t header;
  int ret;

  /* check arguments */
  
  if ( strlen(subsystemname) > MAP_NAMELENGTH-1 )
    {
      map_error(MAP_USER_ERROR_ARGUMENT,
		"map_rem_sub: Subsystem name (%s) longer than %d characters.",
		subsystemname, MAP_NAMELENGTH-1);
      return MAP_USER_ERROR_ARGUMENT;
    }
  
  /* open file for read and write */
  
  fd = map_open_rw(filename);
  if ( fd < MAP_OK )
    {
      map_syserror(fd,
		   "map_rem_sub: can't open (for read/write) TheMap (%s).",
		   filename);
      return fd;
    }
  
  /* find location of subsystem in list */
  
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
      ret = map_find_same(fd,ROOT,subsystemname,0,NULL_LOC,
				    &tindex,&succ_loc,&subsystem_loc);
      if ( ret < 0 )
	{
	  if ( ret == MAP_USER_WARN_NOMATCH )
	    map_warn(ret,
		     "map_rem_sub: subsystem (%s) does not exist",
		     subsystemname);
	  map_close_map(fd);
	  return ret;
	}
      ret = map_find_pred(fd,ROOT,subsystemname,0,NULL_LOC,
			       &tindex,&succ_loc,&pred_loc);
      if ( ret < 0 )
	{
	  map_close_map(fd);
	  return ret;
	}
      status = map_read_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
      
      /* unlink entry */
      
      /* update pointer to new entry */
      
      status = map_overwrite_pointer(&subsystem.next,sizeof subsystem.next,fd,pred_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
    }
  else

    /* version 2 */

    {
      ret = map_find_same(fd,NULL_LOC,subsystemname,
			  header.table_used,
			  header.table,&tindex,&succ_loc,
				    &subsystem_loc);
      if ( ret < 0 )
	{
	  if ( ret == MAP_USER_WARN_NOMATCH )
	    map_warn(ret,
		     "map_rem_sub: subsystem (%s) does not exist",
		     subsystemname);
	  map_close_map(fd);
	  return ret;
	}
      if ( header.table_used-tindex > 1 )
	{
	  /* read portion of table after entry and rewrite it over entry */

	  table = (nametable_t *) calloc(header.table_used-tindex,sizeof(nametable_t));
	  status = map_read_nametable(table,(header.table_used-tindex)*sizeof(nametable_t),fd,header.table+(tindex+1)*(int)sizeof(nametable_t));
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	  status = map_overwrite_nametable(table,(header.table_used-tindex)*sizeof(nametable_t),fd,header.table+tindex*(int)sizeof(nametable_t));
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	  free(table);
	}
      header.table_used--;
      status = map_overwrite_subtableheader(&header,sizeof header,fd,header_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
    }

  /* close the file */
  
  return map_close_map(fd);
  
}
