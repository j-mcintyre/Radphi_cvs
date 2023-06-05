/* 
        map_overwrite_char.c
         Created           :  12-Jan-2000 by Richard Jones
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <stdlib.h>
#include <unistd.h>
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
 
int map_overwrite_char(const char *s, size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
    char *my_s=(char *)s;

#ifdef LITTLE_ENDIAN
#endif

    retval = map_overwrite(my_s, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
#endif

    return retval;
 
}
