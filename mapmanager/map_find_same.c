/* 
   map_find_same.c
   Created           : 16-OCT-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                      */
#include <stdlib.h>
#include <string.h>
#include <map_manager.h> 
#include <map_internal.h>

/***********************************************************************/
/*                                                                     */
/*   MAP_FIND_SAME                                                     */
/*   -------------                                                     */
/*                                                                     */
/*         Created     : 16-OCT-1992    Author : Thom Sulanke          */
/*         Purpose     : find location of name in file.  If name is    */
/*             		 not in file, return MAP_USER_WARN_NOMATCH     */
/*                                                                     */
/***********************************************************************/
 
int map_find_same(int fd, pointer_t pointer, const char name[], 
		  int table_used, 
		  pointer_t table_loc, int *tindex, pointer_t *succ_loc,
		  pointer_t *adr)
     
{
  pointer_t same_loc;
  name_t entry;
  int status;
  nametable_t *table;
  
  if ( table_loc == NULL_LOC )

    /* version 1 */
    
    {
      *tindex = -2;
      *succ_loc = NULL_LOC;
      status = map_read_pointer(&same_loc,sizeof same_loc,fd,pointer);
      if ( status != MAP_OK ) return status;
      while ( same_loc != NULL_LOC )
	{
	  status = map_read_name(&entry,sizeof entry,fd,same_loc);
	  if ( status != MAP_OK ) return status;
	  *succ_loc = entry.next;
	  if ( strncmp(name,entry.name,MAP_NAMELENGTH) == 0 )
	    break;
	  same_loc = entry.next;
	}  /* While */
      
    }
  else

    /* version 2 */
    
    { 
      *tindex = -1;
      *succ_loc = NULL_LOC;
      same_loc  = NULL_LOC;
      if ( table_used > 0 )
 	{
	  table = (nametable_t *) calloc(table_used,sizeof(nametable_t));
	  status = map_read_nametable(table,table_used*sizeof(nametable_t),fd,table_loc);
	  if ( status != MAP_OK )
	    {
	      free(table);
	      map_close_map(fd);
	      return status;
	    }
	  *tindex = 0;
	  while ( *tindex < table_used && strncmp(name,table[*tindex].name,MAP_NAMELENGTH) != 0 )
	    (*tindex)++;
	  if ( *tindex < table_used )
	    same_loc = table[*tindex].loc;

	  if ( *tindex+1 < table_used )
	    *succ_loc = table[*tindex+1].loc;

	  if ( *tindex == table_used )
	    *tindex = -1;

	  free(table);
	}
      
    }

  if ( same_loc == NULL_LOC )
    return MAP_USER_WARN_NOMATCH;
  else{
    *adr = same_loc;
    return MAP_OK;
  }

}






