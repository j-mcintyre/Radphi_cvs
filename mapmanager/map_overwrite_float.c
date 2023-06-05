/* 
        map_overwrite_float.c
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
 
int map_overwrite_float(const float *f, size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
    float *my_f;
 
#ifdef LITTLE_ENDIAN
    my_f = (float *)malloc(nbytes);
    memcpy(my_f,f,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_f,nbytes/sizeof(uint32));
#else
    my_f = f;
#endif

    retval = map_overwrite(my_f, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    free(my_f);
#endif

    return retval;
 
}
