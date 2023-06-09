/* 
        map_write.c
         Created           :  9-OCT-1992 by Thom Sulanke
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
/*         Created     :  9-OCT-1992    Author : Thom Sulanke          */
/*         Purpose     : add to end of file                            */
/*                                                                     */
/***********************************************************************/
 
int   map_write(const void *chunck, size_t nbytes, int fd, pointer_t *adr)
 
{
 
  pointer_t end_before_write;

/* set file position to end of file */
  
  if ( (end_before_write = (int)lseek(fd, 0, SEEK_END)) == -1 )
    {
      map_syserror(MAP_SYSTEM_ERROR_IO,
		   "error positioning to write to end of TheMap");
      map_close_map(fd);
      return MAP_SYSTEM_ERROR_IO;
    }
  
  /* do actual write */
  
  if ( write(fd, chunck, nbytes) == -1 )
    {
      map_syserror(MAP_SYSTEM_ERROR_IO,
		   "error writing to end of TheMap");
      map_close_map(fd);
      return MAP_SYSTEM_ERROR_IO;
    }
  
  *adr = end_before_write;
  return MAP_OK;
    
}
