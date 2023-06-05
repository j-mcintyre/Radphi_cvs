/*
 * -R.T. Jones  July 30, 2003
 *  Added support for ntuple chains longer than 2^31 events.
 *  To do that we had to change a number of variables and arguments
 *  from int to long long int, and hoped not to break anything.  So
 *  far we believe that this should work on piaf without any changes
 *  to paw itself.  To disable, comment out the following line.
 */
#define LONG_LONG_CHAIN 1
/*
 *  pawchn.h  --
 *	Map the /PAWCHN/, /PAWCHC/, /PAWCH2/ and /PAWCH3/ commons
 *
 *  Original:  3-Apr-1995 11:56
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: c_pawchn.h,v 1.1.1.1 2004/03/20 20:05:12 radphi Exp $
 *
 *  $Log: c_pawchn.h,v $
 *  Revision 1.1.1.1  2004/03/20 20:05:12  radphi
 *  cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
 *
 *  Revision 1.3  1999/06/28 15:08:28  couet
 *  - use now cfortran.h in $CVSCOSRC
 *
 *  Revision 1.2  1996/04/23 18:37:49  maartenb
 *  - Add RCS keywords
 *
 *
 */

#ifndef CERN_PAWCHN
#define CERN_PAWCHN

#include	<cfortran/cfortran.h>


typedef struct {
	int	chain;
#ifdef LONG_LONG_CHAIN
	int	padding;
	long long int	nchevt;
	long long int	ichevt;
#else
	int	nchevt;
	int	ichevt;
#endif
} pawchn_def;

#define PAWCHN COMMON_BLOCK(PAWCHN,pawchn)
COMMON_BLOCK_DEF(pawchn_def,PAWCHN);


typedef struct {
	char	cfile[128];
} pawchc_def;

#define PAWCHC COMMON_BLOCK(PAWCHC,pawchc)
COMMON_BLOCK_DEF(pawchc_def,PAWCHC);


typedef struct {
	char	curchn[80];
} pawch2_def;

#define PAWCH2 COMMON_BLOCK(PAWCH2,pawch2)
COMMON_BLOCK_DEF(pawch2_def,PAWCH2);


typedef struct {
	int	chnxmt;
} pawch3_def;

#define PAWCH3 COMMON_BLOCK(PAWCH3,pawch3)
COMMON_BLOCK_DEF(pawch3_def,PAWCH3);


#endif	/*	CERN_PAWCHN	*/
