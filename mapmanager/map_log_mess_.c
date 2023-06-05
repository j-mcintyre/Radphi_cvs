/* 
        map_log_mess_.c
         Created           :  29-MAR-1994 by Thom Sulanke
*/

#include "map_manager.h"
#include "map_internal.h"

/***********************************************************************/
/*                                                                     */
/*   map_log_mess_.c                                                    */
/*   ----------                                                        */
/*                                                                     */
/*         Created     :  23-Mar-1994    Author : Thom Sulanke          */
/*         Purpose     :  interface routine to be called from FORTRAN  */
/*                                                                     */
/***********************************************************************/
 
int map_log_mess_( int *error_code, int *log )
 
{

    return map_log_mess(*error_code,*log);
}
