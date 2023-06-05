/* 
        map_write_char.c
         Created           :  12-Jan-2000 by Richard Jones
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <stdlib.h>
#include <unistd.h>
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
 
int map_write_char(const char *s, size_t nbytes, int fd, pointer_t *adr)
 
{
    int retval;
    char *my_s=(char *)s;

#ifdef LITTLE_ENDIAN
#endif

    retval = map_write(my_s, nbytes, fd, adr);

#ifdef LITTLE_ENDIAN
#endif

    return retval;
 
}
