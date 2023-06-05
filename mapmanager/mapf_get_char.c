/* 
	mapf_get_char.c     : 09-Jul-1994 by R. Hackenburg,
	                    : from Thom Sulanke's mapmanager
	                    : C routines.
*/

#include <stdio.h>                   /* I/O definitions                       */
#include <string.h>
#include "map_manager.h"
#include "map_internal.h"
#define TRUE  1
#define FALSE 0

/*  Prototypes for routines from strbytec_*.c   */

/*                        input        input          output       input      output  */
	void strfc( char *Fstring, int Flength, char *Cstring, int Clength, int *lnb );

/*                        input        input          input        output     output  */
	void strcf( char* Cstring, int Clength, int Flength, char* Fstring, int* lnb );

/***********************************************************************/
/*                                                                     */
/*   mapf_get_char                                                     */
/*                                                                     */
/***********************************************************************/
 
#ifdef vaxc
int mapf_get_char( char *filename, char *subsystemname
                 , char *itemname, char *string
                 , int *atime, int *firsttime, int *stringLength
	         , int filename_len, int subsystemname_len, int itemname_len, int string_len )
#else
int mapf_get_char_( char *filename, char *subsystemname
                  , char *itemname, char *string
                  , int *atime, int *firsttime, int *stringLength
	          , int filename_len, int subsystemname_len, int itemname_len, int string_len )
#endif
/*
	FORTRAN Callable, as:

	LOGICAL FUNCTION mapf_get_char( Filename, Subsystem, Itemname, String, Atime, Firsttime, StringLength )

*  Inputs:
	CHARACTER*(*) Filename       !Map file name.
	CHARACTER*(*) Subsystemname  !Map subsytem name.
	CHARACTER*(*) Itemname       !Map item name.

*  Output:
	CHARACTER*(*) String         !String, to be filled from the Map.

*  Input:
	INTEGER       Atime          !Map "atime" -- time requested.

*  Outputs:
	INTEGER       Firsttime      !Map "firsttime" -- time of items retrieved.
	INTEGER       StringLength   !Length of Array (no. of elements).

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

	status = map_get_char( filename_copy, subsystemname_copy, itemname_copy, string_len
	                     , string, *atime, firsttime, stringLength );

/*	status non-zero means it failed.     */
	if ( status ) { return( FALSE ); }

        return( TRUE );
}

