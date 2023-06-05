/* 
        map_overwrite_name.c
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
 
int map_overwrite_name(const name_t *name, size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
    name_t *my_name;
 
#ifdef LITTLE_ENDIAN
    my_name = (name_t *)malloc(nbytes);
    memcpy(my_name,name,nbytes);
    endian_convertType1_uint32((uint32 *)my_name,nbytes);
#else
    my_name = name;
#endif

    retval = map_overwrite(my_name, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    free(my_name);
#endif

    return retval;
 
}
