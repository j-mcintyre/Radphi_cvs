/* 
        map_add_item.c
         Created           :  10-Nov-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <string.h>
#include "map_manager.h"
#include "map_internal.h"
 

/***********************************************************************/
/*                                                                     */
/*   MAP_ADD_ITEM                                                      */
/*   -----------------                                                 */
/*                                                                     */
/*         Created     :  10-NOV-1992    Author : Thom Sulanke         */
/*         Purpose     :  add a new item to map                        */
/*                                                                     */
/***********************************************************************/
 
int map_add_item( const char filename[], const char subsystemname[], 
            const char itemname[], int length, int type )
{
  /* call internal version with default table size */

  return map_add_item1( filename, subsystemname, 
		       itemname, length, type, INIT_TABLE_SIZE );

}
