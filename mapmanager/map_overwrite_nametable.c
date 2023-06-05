/* 
        map_overwrite_nametable.c
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
/*   MAP_overwrite                                                     */
/*   -------------                                                     */
/*                                                                     */
/*         Created     :  12-Jan-2000   Author : Richard Jones         */
/*         Purpose     : append a structure to the end of the map      */
/*                                                                     */
/***********************************************************************/
 
int map_overwrite_nametable(const nametable_t *nametable, size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
    nametable_t *my_nametable;
    int entry, entries;
 
#ifdef LITTLE_ENDIAN
    my_nametable = (nametable_t *)malloc(nbytes);
    memcpy(my_nametable,nametable,nbytes);
    entries = nbytes/sizeof(nametable_t);
    for (entry=0; entry<entries; entry++) {
        endian_convertType1_uint32((uint32 *)&my_nametable[entry].loc,
                                   sizeof(pointer_t));
    }
#else
    my_nametable = nametable;
#endif

    retval = map_overwrite(my_nametable, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    free(my_nametable);
#endif

    return retval;
 
}
