/* 
        map_write_subtableheader.c
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
 
int map_write_subtableheader(const subtableheader_t *subtableheader, size_t nbytes, int fd, pointer_t *adr)
 
{
    int retval;
    subtableheader_t *my_subtableheader;
 
#ifdef LITTLE_ENDIAN
    my_subtableheader = (subtableheader_t *)malloc(nbytes);
    memcpy(my_subtableheader,subtableheader,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_subtableheader,
                               nbytes/sizeof(uint32));
#else
    my_subtableheader = subtableheader;
#endif

    retval = map_write(my_subtableheader, nbytes, fd, adr);

#ifdef LITTLE_ENDIAN
    free(my_subtableheader);
#endif

    return retval;
 
}
