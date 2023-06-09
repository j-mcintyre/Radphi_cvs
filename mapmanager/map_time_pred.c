/* 
        map_time_pred.c
         Created           : 16-NOV-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                       */
#include <stdlib.h>
#include <map_manager.h>
#include <map_internal.h>


/***********************************************************************/
/*                                                                     */
/*   MAP_TIME_PRED                                                     */
/*   -------------                                                     */
/*                                                                     */
/*         Created     : 16-NOV-1992    Author : Thom Sulanke          */
/*         Purpose     : find predessor (in list, successor in time)   */
/*                          of array of values, also return index of   */
/*                          predessor. An index of -1 indicates no     */
/*                          predessor. An index of -2 indicates no     */
/*                          table.                                     */
/*                                                                     */
/***********************************************************************/
 
int map_time_pred( int fd, pointer_t pointer, item_t item, int atime, 
		  int *tindex, pointer_t *suc_loc, int *suc_time,
		  pointer_t *adr)
 
{
    pointer_t pred_loc, next_loc;
    arrayheader_t entry;
    int status;
    table_t *table;
    
    if ( item.list < 0 )

/* use table */

      {
	if ( item.table_used == 0 )
	  *tindex = -1;
	else
	  {
	    table = (table_t *) calloc(item.table_used,sizeof(table_t));
	    status = map_read_table(table,item.table_used*sizeof(table_t),fd,item.table);
	    if ( status != MAP_OK )
	      {
		free(table);
		map_close_map(fd);
		return status;
	      }
	    *tindex = -1;
	    while ( *tindex < (item.table_used)-1 && table[*tindex+1].time > atime )
	      (*tindex)++;
	  }

	if ( *tindex == -1 )
	  pred_loc = pointer;	
	else
	  pred_loc = table[*tindex].loc;
	
	if ( *tindex == item.table_used-1 )
	  {
	    *suc_loc = NULL_LOC;
	    *suc_time = 0;
	  }
	else
	  {
	    *suc_loc = table[*tindex+1].loc;
	    *suc_time = table[*tindex+1].time;
	  }
	if ( item.table_used != 0 )
	  free(table);
      }

    else

/* use linked list */

      {
	*tindex = -2;
	pred_loc = pointer;
	next_loc = item.list;
	
	while ( next_loc != NULL_LOC )
	  {
	    status = map_read_arrayheader(&entry,sizeof entry,fd,next_loc);
	    if ( status != MAP_OK )
	      {
		map_close_map(fd);
		return status;
	      }
	    if ( entry.time <= atime )
	      break;
	    pred_loc = next_loc;
	    next_loc = entry.next;
	  }  /* While */
	*suc_loc = next_loc;
	*suc_time = entry.time;
      }
   
    *adr = pred_loc;
    return MAP_OK;
    
  }
