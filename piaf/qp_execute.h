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
 *  qp_execute.h  --
 *	Declarations for the query interpreter, and
 *	data structures.
 *
 *  Original: 23-Nov-1994 17:06
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: qp_execute.h,v 1.1.1.1 2004/03/20 20:05:12 radphi Exp $
 *
 *  $Log: qp_execute.h,v $
 *  Revision 1.1.1.1  2004/03/20 20:05:12  radphi
 *  cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
 *
 *  Revision 1.4  1996/04/23 18:38:38  maartenb
 *  - Add RCS keywords
 *
 *
 */

#ifndef CERN_QP_EXECUTE
#define CERN_QP_EXECUTE


#include	"qp_query.h"
#include	"qp_command.h"

void
qp_execute(
	QueryExe *	qxp,
#ifdef LONG_LONG_CHAIN
	long long int	start_evt,
	long long int	nevt,
#else
	long		start_evt,
	long		nevt,
#endif
	QPCmd		*cmd,
	int *		errp
	);

#endif	/*	CERN_QP_EXECUTE	*/
