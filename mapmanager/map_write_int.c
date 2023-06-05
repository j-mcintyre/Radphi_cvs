/* 
        map_write_int.c
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
 
int map_write_int(const int *i, size_t nbytes, int fd, pointer_t *adr)
 
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

    retval = map_write(my_i, nbytes, fd, adr);

#ifdef LITTLE_ENDIAN
    free(my_i);
#endif

    return retval;
 
}
