/* 
        map_write_arrayheader.c
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
 
int map_write_arrayheader(const arrayheader_t *arrayheader,
                          size_t nbytes, int fd, pointer_t *adr)
 
{
    int retval;
    arrayheader_t *my_arrayheader;
 
#ifdef LITTLE_ENDIAN
    my_arrayheader = (arrayheader_t *)malloc(nbytes);
    memcpy(my_arrayheader,arrayheader,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_arrayheader,nbytes/sizeof(uint32));
#else
    my_arrayheader = arrayheader;
#endif

    retval = map_write(my_arrayheader, nbytes, fd, adr);

#ifdef LITTLE_ENDIAN
    free(my_arrayheader);
#endif

    return retval;
 
}
