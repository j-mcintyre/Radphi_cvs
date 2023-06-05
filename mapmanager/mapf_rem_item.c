/* 
	mapf_rem_item.c    : 09-Feb-1994 by R. Hackenburg,
	                   : from Thom Sulanke's mapmanager
	                   : C routines.
*/

#include <stdio.h>                   /* I/O definitions                       */
#include <string.h>
#include "map_manager.h"
#include "map_internal.h"
#define TRUE  1
#define FALSE 0

/*  Prototype for a routine from strbytec_vms.c   */

	void strfc( char *filename, int filename_len, char *filename_copy, int filename_copy_len, int *filename_lnb );

/***********************************************************************/
/*                                                                     */
/*   mapf_add_item                                                     */
/*                                                                     */
/***********************************************************************/


#ifdef vaxc
int mapf_rem_item(  char *filename, char *subsystemname, char *itemname
	          , int filename_len, int subsystemname_len, int itemname_len )
#else
int mapf_rem_item_( char *filename, char *subsystemname, char *itemname
	          , int filename_len, int subsystemname_len, int itemname_len )
#endif
/*
	FORTRAN Callable, as:

	LOGICAL FUNCTION mapf_rem_item( Filename, Subsystem, Itemname )

*  Inputs:
	CHARACTER*(*) Filename       !Map file name.
	CHARACTER*(*) Subsystemname  !Map subsytem name.
	CHARACTER*(*) Itemname       !Map item name.

*  Returns:
*	.TRUE. for success,
*	.FALSE. for general map failure.

*  Description:
	This is an interface for SGI FORTRAN to Thom Sulanke's
	mapmanager C routines, and consists of somewhat modified
	versions of those originally written by Thom Sulanke.

	Note that the names are changed here, to resolve an inconsistancy
	between DEC and SGI with respect to FORTRAN global naming conventions.
	An appended "_", which is implicit (and hidden) under SGI, is
	present here, since SGI does it that way.  To prevent the
	mapmanager routines and their respective interfaces from having
	the same names under VMS, where DEC doesn't do it that way (ie, no
	hidden, appended "_"), "map_" has been changed to "mapf_" in all cases,
	DEC and SGI alike.

*/
 
{
    char filename_copy[MAP_FILENAMELENGTH];
    char subsystemname_copy[MAP_NAMELENGTH];
    char itemname_copy[MAP_NAMELENGTH];
    int  filename_lnb;
    int  subsystemname_lnb;
    int  itemname_lnb;

    int  status;

/*	Copy the FORTRAN-defined filename into a local, C-defined one:  */
	strfc( filename, filename_len, filename_copy, MAP_FILENAMELENGTH, &filename_lnb );

/*	Copy the FORTRAN-defined itemname into a local, C-defined one:  */
	strfc( itemname, itemname_len, itemname_copy, MAP_NAMELENGTH, &itemname_lnb );

/*	Copy the FORTRAN-defined subsystemname into a local, C-defined one:  */
	strfc( subsystemname, subsystemname_len, subsystemname_copy, MAP_NAMELENGTH, &subsystemname_lnb );

	status = map_rem_item(filename_copy,subsystemname_copy,itemname_copy);

/*	status non-zero means it failed.     */
	if ( status ) { return( FALSE ); }
        else          { return( TRUE  ); }

}
