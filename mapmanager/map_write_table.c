/* 
        map_write_table.c
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
 
int map_write_table(const table_t *table, size_t nbytes, int fd, pointer_t *adr)
 
{
    int retval;
    table_t *my_table;
 
#ifdef LITTLE_ENDIAN
    my_table = (table_t *)malloc(nbytes);
    memcpy(my_table,table,nbytes);
    endian_convertTypeN_uint32((uint32 *)my_table,nbytes/sizeof(uint32));
#else
    my_table = table;
#endif

    retval = map_write(my_table, nbytes, fd, adr);

#ifdef LITTLE_ENDIAN
    free(my_table);
#endif

    return retval;
 
}
