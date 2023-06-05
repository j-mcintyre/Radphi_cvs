/*
 *  cern_types.h  --
 *	Separate out the type declaration for storage
 *	size sensitive quantities.
 *
 *  Original: 25-Nov-1994 09:53
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: cern_types.h,v 1.1 2004/07/01 17:59:20 radphi Exp $
 *
 *  $Log: cern_types.h,v $
 *  Revision 1.1  2004/07/01 17:59:20  radphi
 *  Makefile
 *      - changed -m486 switch to -march=i486 to make gcc 3.3 happy [rtj]
 *  buildtree.tar
 *      - added a link to piafs in build tree, for clarity [rtj]
 *  npantup.c, qp_execute.c
 *      - fixed a divide-by-zero error in the case where PCNTWK.nchain==0
 *        and the progress dialog attempts to use it as an increment [rtj]
 *  qp_hbook_if.h, piafront.c, cern_types.h, smap.h, piafc.c
 *      - add #ifndef / #endif protections around defines of constants found
 *        in standard POSIX headers, which are being redefined [rtj]
 *
 *  Revision 1.7  1996/04/23 18:37:54  maartenb
 *  - Add RCS keywords
 *
 *
 */

#ifndef CERN_CERN_TYPES
#define CERN_CERN_TYPES

#include	<float.h>
#include	<limits.h>

#include	"bool.h"


typedef char *			String;

#define	QP_STR_MAX	32

typedef struct _string32_ {
	char	s[QP_STR_MAX];
} String32;


typedef unsigned int		UInt32;

#ifndef UINT32_MAX
#define	UINT32_MAX	UINT_MAX
#endif


/* typedef unsigned long long	UInt64; */
typedef unsigned long 		UInt64;

#ifndef UINT64_MAX
#define	UINT64_MAX	ULONG_MAX
#endif


typedef int			Int32;

#ifndef INT32_MAX
#define	INT32_MAX	INT_MAX
#endif
#ifndef INT32_MIN
#define	INT32_MIN	INT_MIN
#endif


/* typedef long long		Int64; */
typedef long			Int64;

#ifndef INT64_MAX
#define	INT64_MAX	LONG_MAX
#endif
#ifndef INT64_MIN
#define	INT64_MIN	LONG_MIN
#endif


typedef float			Float32;

#ifndef FLOAT32_MAX
#define	FLOAT32_MAX	FLT_MAX
#endif
#ifndef FLOAT32_MIN
#define	FLOAT32_MIN	FLT_MIN
#endif

#ifndef	FLOAT_INT_LIMIT
#define	FLOAT_INT_LIMIT	1.0e10	/* a garanteed minimal value with no fraction */
#endif


typedef double			Float64;

#ifndef FLOAT64_MAX
#define	FLOAT64_MAX	DBL_MAX
#endif
#ifndef FLOAT64_MIN
#define	FLOAT64_MIN	DBL_MIN
#endif

#ifndef DOUBLE_INT_LIMIT
#define	DOUBLE_INT_LIMIT	1.0e20	/* a garanteed minimal value with no fraction */
#endif


#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#endif	/*	CERN_CERN_TYPES	*/
