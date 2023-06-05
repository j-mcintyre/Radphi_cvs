/* 
   map_find_item.c
   Created           : 16-NOV-1992 by Thom Sulanke
*/

 
#include <stdio.h>                   /* I/O definitions                      */
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_FIND_ITEM                                                     */
/*   -------------                                                     */
/*                                                                     */
/*         Created     : 16-NOV-1992    Author : Thom Sulanke          */
/*         Purpose     : return location of item within map file       */
/*                                                                     */
/***********************************************************************/
 
int map_find_item( int fd, const char subsystemname[], 
			const char itemname[],pointer_t *adr )
     
{
  pointer_t subsystem_loc, item_loc, succ_loc, header_loc;
  subtableheader_t header;
  int status;
  subsystem_t subsystem;
  int tindex;
  
  /* find location of subsystem */
  
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
    status = map_find_same(fd,ROOT,subsystemname,0,NULL_LOC,&tindex,
			   &succ_loc,&subsystem_loc);
  else
    status = map_find_same(fd,NULL_LOC,subsystemname,
				  header.table_used,header.table,&tindex,
				  &succ_loc,&subsystem_loc);
  
  if ( status < 0 )
    {
      map_warn(status,
		"no such subsystem %s",
		subsystemname);

      *adr = 0;
      return status; 
    }
  
  /* find location in list of entry */
  
  status = map_read_subsystem(&subsystem,sizeof subsystem,fd,subsystem_loc);
  if ( status != MAP_OK )
    {
      map_close_map(fd);
      return status;
    }
  
  if ( subsystem.list >= 0 )
    
    /* version 1 */
    
    status = map_find_same(fd, subsystem_loc+SUBSYSTEM_LIST_OFFSET, 
			     itemname, 
			     0, NULL_LOC, &tindex, &succ_loc,&item_loc);

    /* version 1 */
    
  else
    status = map_find_same(fd, NULL_LOC, itemname,
			       subsystem.table_used, subsystem.table, &tindex,
			       &succ_loc,&item_loc);

  if ( status < 0 )
    {
      map_warn(status,
		"no such item (%s) in subsystem (%s)",
		itemname,
		subsystemname);

      *adr = 0;
      return status;
    }
  *adr = item_loc;
  return MAP_OK; 
}
