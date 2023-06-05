/* 
        map_stat_item.c
         Created           : by Const Olchanski
*/

#include <stdio.h>                   /* I/O definitions                       */
#include <string.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_STAT_ITEM                                                     */
/*   -------------                                                     */
/*                                                                     */
/*         Created     : 06-NOV-1993    Author : Const Olchanski       */
/*         Purpose     : stat the given item                           */
/*                                                                     */
/***********************************************************************/
 
int map_stat_item( const char filename[], const char subsystemname[], const char itemname[],
		  int *length, int *type)
 
{
    int fd ;
    pointer_t item_loc;
    item_t item;
    int status;

/* set default for early return */

    *length = 0;
    *type = 0;

/* open file for read */
 
    fd = map_open_ro(filename);
    if ( fd < 0 )
    {
        map_syserror(fd,
		     "map_stat_item: can't open Map file (%s)",
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


    *length = item.length;
    *type = item.type;

    return map_close_map(fd);
 
}
