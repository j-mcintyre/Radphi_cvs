/* 
        map_read_nametable.c
         Created           :  12-Jan-2000 by Richard Jones
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <stdlib.h>
#include <unistd.h>
#include <map_manager.h>
#include <map_internal.h>
 

/***********************************************************************/
/*                                                                     */
/*   MAP_READ                                                          */
/*   --------                                                          */
/*                                                                     */
/*         Created     :  12-Jan-2000   Author : Richard Jones         */
/*         Purpose     : read a structure from the map                 */
/*                                                                     */
/***********************************************************************/
 
int map_read_nametable(nametable_t *nametable,
                       size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
    int entry, entries;
 
    retval = map_read(nametable, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    entries = nbytes/sizeof(nametable_t);
    for (entry=0; entry<entries; entry++) {
        endian_convertType1_uint32((uint32 *)&nametable[entry].loc,
                                   sizeof(pointer_t));
    }
#endif

    return retval;
 
}
