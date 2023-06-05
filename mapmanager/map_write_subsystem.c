/* 
        map_write_subsystem.c
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
 
int map_write_subsystem(const subsystem_t *subsystem, size_t nbytes, int fd, pointer_t *adr)
 
{
    int retval;
    subsystem_t *my_subsystem;
 
#ifdef LITTLE_ENDIAN
    my_subsystem = (subsystem_t *)malloc(nbytes);
    memcpy(my_subsystem,subsystem,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_subsystem,nbytes/sizeof(uint32));
    memcpy(my_subsystem->name,subsystem->name,MAP_NAMELENGTH);
#else
    my_subsystem = subsystem;
#endif

    retval = map_write(my_subsystem, nbytes, fd, adr);

#ifdef LITTLE_ENDIAN
    free(my_subsystem);
#endif

    return retval;
 
}
