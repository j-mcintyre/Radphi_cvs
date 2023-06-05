/* 
   map_find_sub.c
   Created           : 16-NOV-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_FIND_SUB                                                      */
/*   ------------                                                      */
/*                                                                     */
/*         Created     : 16-NOV-1992    Author : Thom Sulanke          */
/*         Purpose     : return location of subsystem within map file  */
/*                                                                     */
/***********************************************************************/
 
int map_find_sub( int fd, const char subsystemname[],pointer_t *adr )
     
{
  pointer_t subsystem_loc, succ_loc, header_loc;
  subtableheader_t header;
  int status;
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
    status = map_find_same(fd,ROOT,subsystemname,0,NULL_LOC,&tindex,&succ_loc,&subsystem_loc);
  else
    status = map_find_same(fd,NULL_LOC,subsystemname,
				  header.table_used,header.table,&tindex,&succ_loc,&subsystem_loc);
  
  if ( status < 0 )
    {
      map_warn(status,
		"no such subsystem %s",
		subsystemname);
      return status; 
    }
  
  *adr = subsystem_loc;
  return MAP_OK; 
}
