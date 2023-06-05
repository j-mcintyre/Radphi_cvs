/*
 * -R.T. Jones  July 30, 2003
 *  Added support for ntuple chains longer than 2^31 events.
 *  To do that we had to change a number of variables and arguments
 *  from int to long long int, and hoped not to break anything.  So
 *  far we believe that this should work on piaf without any changes
 *  to paw itself.  To disable, comment out the following line.
 */
#ifndef LONG_LONG_CHAIN
#define LONG_LONG_CHAIN 1
#endif
/* 
 * -R.T. Jones, M. Kornicer  Nov. 7, 2002
 *  A local modification was introduced to prevent piaf slaves from unnecessary
 *  doing useless open/close cycles on files that the slave will never access.  
 *  The ring-buffer semantics of the query processor assume that a consumer of
 *  ntuple data wants to cycle through all of the files in a chain in sequence.
 *  These semantics are implemented in the PCNEXT function, which closes the
 *  current ntuple file and opens the next in the chain.  This is the desired
 *  behavior for local file processing, but for parallel operations it results
 *  in a lot of extra open/close cycles where slaves are simply spacing ahead
 *  to the desired location in the chain, or finding their way back to the
 *  start.  A new function h_next_nt() has been introduced, which slaves can
 *  use as a front-end to PCNEXT, which obeys these semantics but avoids any
 *  actual file operations unless the ntuple rows in the given file are
 *  actually needed.  To do its job, h_next_nt() needs to know how many rows
 *  are stored in each file, which was already discovered, but not stored,
 *  in h_load_nt().  A modification to h_load_nt() has been introduced to
 *  save this information, and in the source of qp_execute.c to use it.
 *  To select these changes, uncomment the following line.
 */
#ifndef PCNEXT_OPTIMIZATION
#define PCNEXT_OPTIMIZATION TRUE
#endif
/*
 *  qp_hbook_if.h  --
 *	Interface of convenience routines for hbook.
 *
 *  Original:  8-May-1995 22:54
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: qp_hbook_if.h,v 1.4 2005/11/01 15:34:30 radphi Exp $
 *
 *  $Log: qp_hbook_if.h,v $
 *  Revision 1.4  2005/11/01 15:34:30  radphi
 *  psexec.F
 *     - just a spelling error in the comments [rtj]
 *  pslvio.F
 *     - make the Format statement more robust [rtj]
 *  qp_execute.c, qp_hbook_if.c qp_hbook_if.h
 *     - fix a bug in qp_exec_init (qp_execute.c) where it was passing a short
 *       row count (PCNEXT_OPTIMIZATION only, in the case where the slave
 *       computes a non-contiguous slice of the ntuple chain) and so reserving
 *       too little space in the ntuple data cache to hold the data [rtj]
 *     - introduced a new function h_nrows_piaf() to return the TOTAL number
 *       of rows in a slave's ntuple processing slice, as opposed to the size
 *       of the first contiguous segment which is returned by h_range_piaf()
 *       and h_next_nt(). [rtj]
 *
 *  Revision 1.3  2005/04/28 19:32:14  radphi
 *  chain.c
 *     - static function find_event() was enclosed in #if 0/#endif, it
 *       is now enabled within the scope of other functions in chain.c [rtj]
 *  qp_execute.c
 *     - fixed a bug in declarations that prevented compilation if LONG_LONG_CHAIN
 *       was not enabled [rtj]
 *     - pass additional arguments ievt,end_evt by reference to h_next_nt()
 *       so that h_next_nt() can update the current row pointer and last row
 *       pointer, as a part of the PCNEXT_OPTIMIZATION feature. [rtj]
 *  qp_hbook_if.h
 *     - modified argument list for h_next_nt() and added new function
 *       h_range_piaf(), both found in qp_hbook_if.c [rtj]
 *  qp_hbook_if.c
 *     - extensive modifications and additions under compilation option
 *       PCNEXT_OPTIMIZATION, to implement a flexible load distribution scheme
 *       between the piaf slaves in the case of ntuple chain processing (see
 *       comments in the code) [rtj]
 *  npantup.c
 *     - modified calc_event_range() function to implement the multi-range
 *       feature of PCNEXT_OPTIMIZATION [rtj]
 *
 *  Revision 1.2  2004/07/01 17:59:21  radphi
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
 *  Revision 1.1.1.1  2004/03/20 20:05:12  radphi
 *  cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
 *
 *  Revision 1.9  1996/09/12 09:28:57  couet
 *  - after NTUPLE commands using a full path name ie:
 *
 *    NT/PLOT //lun1/10.x
 *
 *    The current working directory was changed to //lun1
 *
 *  Revision 1.8  1996/04/23 18:38:44  maartenb
 *  - Add RCS keywords
 *
 *
 */

#ifndef CERN_QP_HBOOK_IF
#define CERN_QP_HBOOK_IF

#include	"cern_types.h"
#include	"hbook_defs.h"
#include	"smap.h"


extern char	h_rwn_chtitl[128];
extern int	h_rwn_nvar;
extern char	h_rwn_tags[MAX_RWN_COLS][9];
extern float	h_rwn_rlow[MAX_RWN_COLS];
extern float	h_rwn_rhigh[MAX_RWN_COLS];


int
h_load_nt(
	char *		id_string,
	char **		id_path,
	int *		idp
	);

#ifdef PCNEXT_OPTIMIZATION
/* new function h_next_nt() defined in qp_hbook_if.c
 * as a kind of replacement for pcnext() in pawlib
 * -rtj-
 */
int
h_next_nt(
	int		id,
# ifdef LONG_LONG_CHAIN
	long long int *	nchrow,
# else
	int *		nchrow,
# endif
	int *		ndim,
	int *		nrow,
	int *		iend,
# ifdef LONG_LONG_CHAIN
	long long int *	next_evt,
	long long int *	end_evt
# else
	int 	      *	next_evt,
	int 	      *	end_evt
# endif
	);

/* new function h_range_piaf() defined in qp_hbook_if.c
 * as a kind of replacement for hrngpf() in packlib/hbook
 * -rtj-
 */
int
h_range_piaf(
	int		slave_count,
	int		slave_index,
# ifdef LONG_LONG_CHAIN
	long long int *	first_evt,
	long long int *	final_evt
# else
	int	      *	first_evt,
	int	      *	final_evt
# endif
	);

int h_nrows_piaf();

#endif

void
h_reset_dir(
	   );

int
h_load_histo(
	char *		id_string,
	int *		idp,
	int *		id_dim
	);

void
h_hnocol(
	int		id,
	int *		colp
	);

#ifdef LONG_LONG_CHAIN
long long int
#else
void
#endif
h_hnoent(
	int		idn,
	bool		use_chain,
	int *		ep
	);

bool
h_flag_1d(
	int		id
	);

bool
h_flag_2d(
	int		id
	);

bool
h_flag_profile(
	int		id
	);


int
h_rwn_getInfo(
	int		idn
);


int
h_rwn_getIndex(
	int		idn,
	char *		const name
);


SMap
h_get_labels(
	int		id,
	char		*chopt
);


void
h_hbook1_labels(
	int		idh,
	char *		title,
	SMap		label_list
);


void
h_hbook2_labels(
	int		idh,
	char *		title,
	SMap		labelx_list,
	SMap		labely_list,
	int		nbin[],
	float		min[],
	float		max[]
);

#endif	/*	CERN_QP_HBOOK_IF	*/
