/* 
   map_dump_.c
   Created           : 29-NOV-1993 by Thom Sulanke
*/

#include <map_manager.h>
#include <map_internal.h>
#include <string.h>

/***********************************************************************/
/*                                                                     */
/*   map_dump_.c                                                       */
/*   -----------                                                       */
/*                                                                     */
/*         Created     :  29-Nov-1993    Author : Thom Sulanke         */
/*         Purpose     :  interface routine to be called from FORTRAN  */
/*                                                                     */
/***********************************************************************/
 
int map_dump_(char *filename, int filename_len)
 
{
    char filename_copy[MAP_FILENAMELENGTH];

    if ( filename_len > MAP_FILENAMELENGTH-1 )
      {
	map_error(MAP_USER_ERROR_ARGUMENT,
		  "MAP_DUMP: filename longer (%d) than %d characters.",
		  filename_len,
		  MAP_FILENAMELENGTH-1);
	return MAP_USER_ERROR_ARGUMENT;
      }
    strncpy(filename_copy,filename,MAP_FILENAMELENGTH);

    return map_dump(filename_copy);
}
