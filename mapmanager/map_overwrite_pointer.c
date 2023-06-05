/* 
        map_overwrite_pointer.c
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
 
int map_overwrite_pointer(const pointer_t *pointer, size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
    pointer_t *my_pointer;
 
#ifdef LITTLE_ENDIAN
    my_pointer = (pointer_t *)malloc(nbytes);
    memcpy(my_pointer,pointer,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_pointer,nbytes/sizeof(uint32));
#else
    my_pointer = pointer;
#endif

    retval = map_overwrite(my_pointer, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    free(my_pointer);
#endif

    return retval;
 
}
