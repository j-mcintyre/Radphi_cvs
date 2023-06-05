/* 
        map_read_subsystem.c
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
 
int map_read_subsystem(subsystem_t *subsystem,
                       size_t nbytes, int fd, pointer_t offset)
 
{
    int retval;
 
    retval = map_read(subsystem, nbytes, fd, offset);

#ifdef LITTLE_ENDIAN
    endian_convertTypeN_uint32((uint32 *)subsystem, nbytes/sizeof(uint32));
    endian_convertTypeN_uint32((uint32 *)subsystem->name,
                               MAP_NAMELENGTH/sizeof(uint32));
#endif

    return retval;
 
}
