/* 
        map_overwrite_int.c
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
 
int map_overwrite_int(const int *i, size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
    int *my_i;
 
#ifdef LITTLE_ENDIAN
    my_i = malloc(nbytes);
    memcpy(my_i,i,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_i,nbytes/sizeof(uint32));
#else
    my_i = i;
#endif

    retval = map_overwrite(my_i, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    free(my_i);
#endif

    return retval;
 
}
