/* 
        map_write_item.c
         Created           :  12-Jan-2000 by Richard Jones
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_WRITE                                                         */
/*   ---------                                                         */
/*                                                                     */
/*         Created     :  12-Jan-2000   Author : Richard Jones         */
/*         Purpose     : append a structure to the end of the map      */
/*                                                                     */
/***********************************************************************/
 
int map_write_item(const item_t *item, size_t nbytes, int fd, pointer_t *adr)
 
{
    int retval;
    item_t *my_item;
 
#ifdef LITTLE_ENDIAN
    my_item = (item_t *)malloc(nbytes);
    memcpy(my_item,item,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_item,nbytes/sizeof(uint32));
    memcpy(my_item->name,item->name,MAP_NAMELENGTH);
#else
    my_item = item;
#endif

    retval = map_write(my_item, nbytes, fd, adr);

#ifdef LITTLE_ENDIAN
    free(my_item);
#endif

    return retval;
 
}
