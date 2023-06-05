/* 
   map_add_sub.c
   Created           :  1-OCT-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_ADD_SUB                                                       */
/*   -----------                                                       */
/*                                                                     */
/*         Created     :  1-OCT-1992    Author : Thom Sulanke          */
/*         Purpose     :  add a new subsystem to map                   */
/*                                                                     */
/***********************************************************************/

int map_add_sub( int fd, const char subsystemname[], int nitems )
{
  pointer_t pred_loc, succ_loc, subsystem_loc, header_loc;
  subsystem_t subsystem, successor;
  int status, tindex;
  nametable_t *table;
  subtableheader_t header;
  
  /* find location in list for entry */
  
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
    status = map_find_pred(fd,ROOT,subsystemname,0,NULL_LOC,&tindex,&succ_loc,&pred_loc);
  else
    /* version 2 */
    status = map_find_pred(fd,NULL_LOC,subsystemname,
			     header.table_used,header.table,&tindex,&succ_loc,&pred_loc);
  
  if ( status < 0 )
    {
      map_close_map(fd);
      return status;
    }
  
  if ( succ_loc != NULL_LOC )
    {
      status = map_read_subsystem(&successor,sizeof successor,fd,succ_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
      if ( strncmp(subsystemname, successor.name, MAP_NAMELENGTH) == 0 )
        {
	  return MAP_OK;
        }
    }
  
  /* fill entry name */
  
  map_fillname(subsystem.name,subsystemname);
  subsystem.list = -VERSION;
  subsystem.table_length = nitems;
  subsystem.table_used = 0;
  table = (nametable_t *) calloc(subsystem.table_length,sizeof(nametable_t));
  status = map_write_nametable(table,subsystem.table_length*sizeof(nametable_t),
			       fd,&subsystem.table);
  free(table);
  if ( status < 0 )
    {
      map_close_map(fd);
      return status;
    }
  
  /* fill in entry next pointer */
  
  if ( header.negversion >= 0 )
    subsystem.next = succ_loc;
  else
    subsystem.next = -VERSION;

  /* insert entry */
  
  status = map_write_subsystem(&subsystem,sizeof subsystem,fd,&subsystem_loc);
  if ( status < 0 )
    return status;
  
  if ( header.negversion >= 0 )
    
    /* version 1 */
    /* update pointer to new entry */
    
    return map_overwrite_pointer(&subsystem_loc,sizeof subsystem_loc,fd,pred_loc);
  
  else
    
    /* version 2 */
    /* add name to table */
    
    {
      if ( header.table_used < header.table_length )
	table = (nametable_t *) calloc(header.table_used+1,sizeof(nametable_t));
      else
	table = (nametable_t *) calloc(header.table_length * TABLE_FACTOR,sizeof(nametable_t));

      if ( header.table_used > 0 )
	{
	  status = map_read_nametable(table,header.table_used*sizeof(nametable_t),
			              fd,header.table);
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	}
      if ( header.table_used-tindex-1 > 0 )
	memmove(&table[tindex+2],&table[tindex+1],
		(header.table_used-tindex-1)*sizeof(nametable_t));
      memcpy(table[tindex+1].name,subsystem.name,MAP_NAMELENGTH);
      table[tindex+1].loc = subsystem_loc;
      header.table_used++;
      if ( header.table_used <= header.table_length )
	{
	  status = map_overwrite_nametable(table,header.table_used*sizeof(nametable_t),
				           fd,header.table);
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	}
      else
	{
	  header.table_length = header.table_length * TABLE_FACTOR;
	  status = map_write_nametable(table,
				       header.table_length*sizeof(nametable_t),fd,&header.table);
	  if ( status < 0 )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	}
      free(table);
      status = map_overwrite_subtableheader(&header,sizeof header,fd,header_loc);
      if ( status != MAP_OK )
	{
	  map_close_map(fd);
	  return status;
	}
    }
  return MAP_OK;
}
