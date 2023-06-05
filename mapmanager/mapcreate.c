/* 
        mapcreate
         Created           : 18-DEC-1992 by Thom Sulanke
*/
 
#include <stdio.h>                   /* I/O definitions                */
#include <map_manager.h>
#include <stdlib.h>

static char  rcsid[] = "$Id: mapcreate.c,v 1.3 2000/01/18 04:15:31 radphi Exp $";
/*
$Log: mapcreate.c,v $
Revision 1.3  2000/01/18 04:15:31  radphi
Changes made by jonesrt@golem
major upgrade to mapamanger package to support little-endian processors;
data on disk remains in big-endian format as before, and is translated on
input and output for little-endian machines (such as the Intel x86).
-rtj

Revision 1.2  1996/11/08 16:21:09  lfcrob
fixed $Log: mapcreate.c,v $
fixed Revision 1.3  2000/01/18 04:15:31  radphi
fixed Changes made by jonesrt@golem
fixed major upgrade to mapamanger package to support little-endian processors;
fixed data on disk remains in big-endian format as before, and is translated on
fixed input and output for little-endian machines (such as the Intel x86).
fixed -rtj
fixed ???

*/



char *progname;


/***********************************************************************/
/*                                                                     */
/*   mapcreate                                                         */
/*   ---------                                                         */
/*                                                                     */
/*         Created     : 18-DEC-1992    Author : Thom Sulanke          */
/*         Purpose     : create a map                                  */
/*                                                                     */
/***********************************************************************/
 
int main( argc, argv )
    int argc ;
    char *argv[]; 
{
    int f1, f2, n ;
    
    progname = argv[0];
    
    if ( argc != 2 )
    {
        printf("Usage: %s themap\n", progname);
	exit(0);
    }
    map_create(argv[1]);

    exit(0);
}
