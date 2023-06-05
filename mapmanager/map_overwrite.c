/* 
        map_overwrite.c
         Created           :  9-OCT-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <stdlib.h>
#include <unistd.h>
#include <map_manager.h>
#include <map_internal.h>

/***********************************************************************/
/*                                                                     */
/*   MAP_OVERWRITE                                                     */
/*   ---------                                                         */
/*                                                                     */
/*         Created     :  9-OCT-1992    Author : Thom Sulanke          */
/*         Purpose     : overwrite portion of file                     */
/*                                                                     */
/***********************************************************************/
 
int map_overwrite(const void *chunck, size_t nbytes, int fd, pointer_t offset)
 
{
 
/* set file position to offset */
 
    if ( lseek(fd, offset, SEEK_SET) == -1 )
    {
        map_syserror(MAP_SYSTEM_ERROR_IO,
		     "error positioning to overwrite in TheMap");
        map_close_map(fd);
        return MAP_SYSTEM_ERROR_IO;
    }

/* do actual write */
 
    if ( write(fd, chunck, nbytes) != nbytes )
    {
        map_syserror(MAP_SYSTEM_ERROR_IO,
		     "error overwriting in TheMap");
        map_close_map(fd);
        return MAP_SYSTEM_ERROR_IO;
    }

    return MAP_OK;
 
}
