/* 
        map_read_pointer.c
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
 
int map_read_pointer(pointer_t *pointer,
                     size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
 
    retval = map_read(pointer, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    endian_convertType1_uint32((uint32 *)pointer, sizeof(pointer_t));
#endif

    return retval;
 
}
