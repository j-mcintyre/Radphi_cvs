/*
 * -R.T. Jones  March 11, 2004
 *  Added support for a user-defined weights that determine the partition
 *  of the ntuple or chain among the slaves.  Default weights of unity are
 *  assigned to each slave at piaf startup.  The overall normalization of
 *  the weights is arbitrary.  The command to assign new weights is:
 *       paw> piaf/message reload <s> <w>
 *  where <s> is a slave index 1..NSLAVE and <w> is an integer weight that
 *  represents the relative number of rows to be processed by that slave.
 *  To disable, comment out the next line.
 */
#define WEIGHTED_NTUPLE_PARTITIONING 1
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
#define PCNEXT_OPTIMIZATION 1
#define DEBUGGING_PCNEXT_OPTIMIZATION 0
/*
 * -R.T. Jones, Nov. 7, 2002
 *  Added a couple of fixes to make things work better under Linux.
 *  To disable these changes, comment out the following line.
 */
#define MAKE_GCC_HAPPY 1
/*
 *  qp_hbook_if.c  --
 *	Implementation of convenience routines for hbook.
 *
 *  Original:  8-May-1995 22:55
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: qp_hbook_if.c,v 1.7 2008/05/29 21:31:00 radphi Exp $
 *
 *  $Log: qp_hbook_if.c,v $
 *  Revision 1.7  2008/05/29 21:31:00  radphi
 *  pfpaw.F, qp_hbook_if.c [rtj]
 *     - disabled the DEBUGGING_*** switches for production running
 *
 *  Revision 1.6  2006/03/30 02:22:04  radphi
 *  npantup.c, qp_execute.c
 *     - changed back to report statistics without rescaling -rtj
 *  qp_hbook_if.c
 *     - bug fix to handle case with only one event per file -rtj
 *     - bug fix to handle switching between ntuples within a given chain -rtj
 *     - other minor fixes -rtj
 *  pfpaw.F
 *     - fixed handling of work dialog boxes to prevent repeated popping -rtj
 *
 *  Revision 1.5  2005/12/30 03:09:08  radphi
 *  qp_hbook_if.c
 *     - turned off DEBUGGING_PCNEXT_OPTIMIZATION by default [rtj]
 *
 *  Revision 1.4  2005/12/30 03:06:47  radphi
 *  qp_hbook_if.c
 *     - moved reset chain_table_len=0 into h_range_piaf() so that changing
 *       chains will start from a fresh table instead of appending [rtj]
 *     - added new debugging prints for DEBUGGING_PCNEXT_OPTIMIZATION [rtj]
 *
 *  Revision 1.3  2005/11/01 15:34:30  radphi
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
 *  Revision 1.2  2005/04/28 19:32:14  radphi
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
 *  Revision 1.1.1.1  2004/03/20 20:05:12  radphi
 *  cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
 *
 *  Revision 1.34  1999/07/05 15:43:35  couet
 *  - hbook_interface.h in now replaced by hbook.h in CVSCOSRC
 *
 *  Revision 1.33  1999/07/01 14:47:01  couet
 *  - cspack_interface.h is replaced by cspack.h in CVSCOSRC
 *
 *  Revision 1.32  1999/06/28 15:09:02  couet
 *  - use now cfortran.h in $CVSCOSRC
 *
 *  Revision 1.31  1997/03/11 11:23:23  couet
 *  - In h_load_nt, the str_del(id_path) are removed because it is done in hntld_C
 *  - In hntld_C the pointer to the string (idpath) is intialized with 0 and the
 *    str_del is done only if the pointer is non zero.
 *
 *    This problems apeared on OSF when nt/print was done on a ntuple which doesn't
 *    exist.
 *
 *  Revision 1.30  1997/03/10 17:03:31  couet
 *  - in h_load_nt iquest(1) was not set at the right place when the ntuple
 *  doesn't exit
 *
 *  Revision 1.29  1996/09/12 09:28:56  couet
 *  - after NTUPLE commands using a full path name ie:
 *
 *    NT/PLOT //lun1/10.x
 *
 *    The current working directory was changed to //lun1
 *
 *  Revision 1.28  1996/09/11 15:01:09  couet
 *  - OLD_PATH is reset in h_load_nt when it exist in errors
 *  - New routine HRSTD to reset OLD_PATH from f77
 *
 *  Revision 1.27  1996/06/27 15:41:32  couet
 *  - Dummy command F/UNIT suppressed (ambiguous with K/UNITS)
 *
 *  Revision 1.26  1996/06/25 15:16:23  couet
 *  - C comment in CVS comments are not allowed !!!
 *
 *  Revision 1.25  1996/06/24 08:50:38  couet
 *  - Bug fix from Dong Zhao <dongzhao@d0sgi0.fnal.gov>:
 *
 *  I got two problems when I used the new Version of PAW/PIAF (96a) from CERN.
 *
 *  First, when I use a large Column Wise Ntuple, (400 variables, 8176 words of
 *  buffer space, 8191 words of record length), it showed:
 *
 *  PAW [1] hi/file 1 /proj5/qcdwz1/zee_1b/v9/data/zee_1b_1214_02.ntp_v9 0
 *  PAW [2] nt/pl 10.z_mass ! 1000
 *   ***** ERROR in HNBUFR : Not enough space in memory : ID=      10
 *
 *  The reason is that in this new version of PAW, HBBUFR reserves space for all
 *  the variables, in my case, 400X(8176+33) > 2000000 PAWC , while in the older
 *  versions, they call a routine called "hgetnt2" in which
 *  "hbname(idn,' ',0,'$CLEAR')" is called to reset the variable addresses.  In the
 *  new version, h_load_nt (corresponds to the old hgetnt) doesn't call hbname to
 *  reset the addresses of the variables.
 *
 *  The fix would be adding a line:
 *
 *          if ( ! PNTOLD.ntold )  HBNAME( idn, " ", 0, "$CLEAR" );
 *
 *  in the end of routine h_load_nt, right before line
 *
 *          * done *
 *
 *  Revision 1.24  1996/05/24 09:14:25  dinofm
 *  Memory leakage on id_path fixed
 *
 *  Revision 1.23  1996/05/17 15:52:12  dinofm
 *  Fixed a bug in h_load_nt. When the ntuple was on Piaf the HGETNT command
 *  was passed to CZPUTA through a common area mapped by cfortran.h.
 *  This caused CZPUTA to crash on AIX, OSF and possibly VMS. I removed any
 *  code related to the mapped common area and just passed the C buffer to
 *  the CZPUTA Fortran wrapper (tested on HP-UX & AIX).
 *
 *  Revision 1.22  1996/05/14 12:23:27  maartenb
 *  - Fix prototypes.
 *
 *  - Fix static bool conversions
 *
 *  Revision 1.21  1996/05/10 12:23:27  dinofm
 *  PIAF debug messages (Ntuple on Piaf etc.) are issued only if the 'verbose'
 *  flag is not 0.
 *
 *  Revision 1.20  1996/05/02 08:06:24  dinofm
 *  Another fix about PIAF chains (ntuple path was used instead than ntuple id
 *  when composing the HGETNT command).
 *
 *  Revision 1.19  1996/04/26 13:21:47  maartenb
 *  - Fixed? the case where we get the number of entries from a chain on piaf
 *
 *  Revision 1.18  1996/04/23 18:38:43  maartenb
 *  - Add RCS keywords
 *
 *
 */

#ifdef MAKE_GCC_HAPPY
# ifdef linux
# define f2cFortran 1
# endif
#endif

#include	<string.h>

#include	"str.h"

#include	<cfortran/cfortran.h>

#include	"c_hcbits.h"
#include	"c_hcbook.h"
#include	"c_hcdire.h"
#include	"c_hcntpar.h"
#include	"c_hcpiaf.h"
#include	"c_pawchn.h"
#include	"c_pawcom.h"
#include	"c_pntold.h"
#include	"c_quest.h"
#include	<cfortran/cspack.h>
#include	"errors.h"
#include	<cfortran/hbook.h>
#include	"kuip_interface.h"
#include	"paw_interface.h"
#include	"qp_hbook_if.h"
#include	"qp_report.h"
#include	"qp_command.h"	/* for MAX_MAP_STRINGS */
#include	"qpflags.h"


#define	CF_TRUE		C2FLOGICAL(1)
#define	CF_FALSE	C2FLOGICAL(0)


char	old_path[MAX_RZ_PATH];

#if defined PCNEXT_OPTIMIZATION
# define MAX_CHAIN_LEN 9999
struct ntuple_chain_file_descriptor_ {
	int	nrows;
	int	slice;
	int	start;
	int	end;
} chain_table[MAX_CHAIN_LEN];
int	chain_table_len;
int	chain_table_id = 0;

  PROTOCCALLSFSUB4(PCHNXT,pchnxt,STRING,INT,STRING,INT)
# define PCHNXT(A1,A2,A3,A4)  CCALLSFSUB4(PCHNXT,pchnxt,STRING,INT,STRING,INT,A1,A2,A3,A4)

# define min(A1,A2)  (((A1)<(A2))? (A1):(A2))
# define max(A1,A2)  (((A1)>(A2))? (A1):(A2))

#endif

int
h_load_nt(
	char *		id_string,
	char **		id_path,
	int *		idp
) {
	char	idtmp_path[MAX_RZ_PATH];
	char	path[MAX_RZ_PATH], *p, *q;
	int	idn, idtmp;
	int	icycle;
	int	ierr;
	int     izero = 0;

#if DEBUGGING_PCNEXT_OPTIMIZATION
	int d;
	FILE *log;
	char filename[99];
	sprintf(filename,"piafs%d.dlog",HCPIAF.mysid);
	log=fopen(filename,"a");
#endif

	/* split string into  path, id and cycle */

	p = str_new( id_string );

	q = strrchr( p, ';' );
	if ( q != 0 ) {
		if ( str_atoi( q+1, 10, &icycle ) ) {
			sf_report( "Illegal cycle '%s'\n", q+1 );
			str_del( p );
			return R_SYNTAX_ERROR;
		}
		*q = '\0';
	} else {
		icycle = 9999;
	}
	q = strrchr( p, '/' );
	if ( q != 0 ) {
		if ( str_atoi( q+1, 10, &idn ) ) {
			sf_report( "Illegal IDN '%s'\n", q+1 );
			str_del( p );
			return R_SYNTAX_ERROR;
		}
		*q = '\0';
		strncpy( path, p, MAX_RZ_PATH-1 );
		path[MAX_RZ_PATH-1] = '\0';
		for( q=path ; *q != '\0' ; q++ ) {
			*q = toupper( *q );
		}
	} else {
		if ( str_atoi( p, 10, &idn ) ) {
			sf_report( "Illegal IDN '%s'\n", p );
			str_del( p );
			return R_SYNTAX_ERROR;
		}
		path[0] = '\0';
	}
	str_del( p );

	PAWCOM.id = idn; /* backward compatibility */
	if ( qp_flags_get( "verbose" ) != 0 ) {
		sf_report( "Loading Ntuple:  path='%s'  id=%d  icycle=%d\n", path,
			idn, icycle );
	}

	/* store the old path */

	HCDIR( old_path, "R" );

	/* cd to new dir / setup chain if necessary */
#if DEBUGGING_PCNEXT_OPTIMIZATION
	fprintf(log,"changing to new dir/chain %s\n",path);
#endif
	PCHNCD( path, ierr );
	if ( ierr == 2 ) {
		QUEST.iquest[0] = -3;
		return R_SYNTAX_ERROR;
	}
	HCDIR( path, " " );
	if ( QUEST.iquest[0] != 0 ) {
		return R_SYNTAX_ERROR;
	}

	/* is the ntuple on piaf ? */
	if ( (HCDIRN.ichtop[HCDIRN.icdir-1] > 1000) &&
		(HCDIRN.ichlun[HCDIRN.icdir-1] != 0) ) {
		HCPIAF.ntpiaf = CF_TRUE;
		if ( qp_flags_get( "verbose" ) != 0 ) {
			sf_report( "Ntuple is on piaf\n" );
		}
	} else {
		HCPIAF.ntpiaf = CF_FALSE;
	}

	/* the id of the ntuple in memory will be : */
	idtmp = idn + PAWCOM.jofset;

	/* Check if the working directory is different from //pawc */
	if ( HCDIRN.ichtop[HCDIRN.icdir-1] != 0 ) {

		/* delete existing copy if necessary */
		if ( HEXIST(idtmp) ) {
			HFIND( idtmp, "h_load_nt" );
			HDCOFL();
			if ( HCBITS.i4 == 0 ) {
				/* idtmp is not an ntuple */
				strcpy( idtmp_path, "???" );
			} else {
				/* get path of idtmp */
				idtmp_path[0] = '\0';
				HGTDIR( idtmp_path );
			}
			if ( 
				( (path[0] != '\0') &&
				(strcmp(path,idtmp_path) != 0) )
			||
				( (path[0] == '\0') &&
				(strcmp(old_path,idtmp_path) != 0) )
			||
				PAWCHN.chain	) {
				HDELET( idtmp );
				HRIN( idn, icycle, PAWCOM.jofset );
			}
		} else {
			/* just load the ntuple into //pawc */
			HRIN( idn, icycle, PAWCOM.jofset );
		}
#if DEBUGGING_PCNEXT_OPTIMIZATION
		fprintf(log,"ntuple %d was just purged!\n",idn);
#endif
	}

	/* setup pointers and flags */
	HFIND( idtmp, "h_load_nt" );
	if ( HCBOOK.lcid == 0 ) {
		/* not found ... */
                HCDIR( old_path, " " );
		QUEST.iquest[0] = -1;
		return R_SYNTAX_ERROR;
	}
	HDCOFL();

	if ( HCBITS.i4 == 0 ) {
		ku_alfa();
		sf_report( "ID %d is not an ntuple.\n", idn );
                HCDIR( old_path, " " );
		QUEST.iquest[0] = -2;
		return R_SYNTAX_ERROR;
	}

	idtmp_path[0] = '\0';
	HGTDIR( idtmp_path );
	*id_path =str_new( idtmp_path );
	*idp = idn;

	HCDIR( path, "R" );	/* retrieve path, because of chain ?? */

	PNTOLD.ntold = ! HNTNEW( idtmp );

	/* check for needed recover */

	if ( PNTOLD.ntold ) {
		if ( PAWC.lq_eqv[HCBOOK.lcid-1-1] <= 0 ) {
			HDELET( idtmp );
			sf_report( "ID %d: Bad header. Try command RECOVER\n",
				idn );
                        HCDIR( old_path, " " );
			QUEST.iquest[0] = -4;
			return R_SYNTAX_ERROR;
		}
	} else {
		if ( PAWC.lq_eqv[HCBOOK.lcid+3-1] <= 0 ) {
			HDELET( idtmp );
			sf_report( "ID %d: Bad header. Try command RECOVER\n",
				idn );
                        HCDIR( old_path, " " );
			QUEST.iquest[0] = -4;
			return R_SYNTAX_ERROR;
		}
	}

	/* determine number of rows (chains/piaf !) */

	if ( PAWCHN.chain && ! HCPIAF.ntpiaf ) {
#ifdef LONG_LONG_CHAIN
		long long int zero = 0;
		int two = 2;
		int len;
		len = strlen(path);
		PAWCHN.nchevt = pchevt_( path, &len, &idn, &zero, &two);
#else
		PAWCHN.nchevt = PCHEVT( path, strlen(path), idn, 0, 2);
#endif

#if PCNEXT_OPTIMIZATION
		if (idn != chain_table_id) {
                        chain_table_id = idn;
                	PAWCHN.nchevt = 0;
                }
#endif

		if ( PAWCHN.nchevt == 0 ) {
#ifdef LONG_LONG_CHAIN
			long long int nchrow;
			int	nrow, iend;
#else
			int	nchrow, nrow, iend;
#endif
#if PCNEXT_OPTIMIZATION
			chain_table_len = 0;
#endif

			nchrow = 0;
			nrow = PAWC.iq_eqv[HCBOOK.lcid+3-1];
			for ( iend=0 ; iend == 0 ; ) {
				int	ncol;
#if defined PCNEXT_OPTIMIZATION
                                int slice;
				chain_table[chain_table_len].nrows = nrow;
				chain_table[chain_table_len].start = 1;
				chain_table[chain_table_len].end = 0;
				char *p = PAWCHC.cfile;
				if (p = strrchr(p,':')) {
				    if (! sscanf(p,":%d",&slice))
					slice = 0;
				    slice = (slice < 0)? 0 : slice;
				    slice = (slice < HCPIAF.ngsize)?
					     slice : HCPIAF.ngsize;
				}
				else {
				    slice = 0;
				}
				chain_table[chain_table_len].slice = slice;
				++chain_table_len;
#endif
#if defined LONG_LONG_CHAIN
				pcnext_( &idtmp, &nchrow, &ncol, &nrow, &iend );
#else
				PCNEXT( idtmp, nchrow, ncol, nrow, iend );
#endif
			}
		}
	}

	if ( PAWCHN.chain && HCPIAF.ntpiaf ) {

#ifdef LONG_LONG_CHAIN
		long long int zero = 0;
		int two = 2;
		int len;
		len = strlen(path);
		PAWCHN.nchevt = pchevt_( path, &len, &idn, &zero, &two);
#else
		PAWCHN.nchevt = PCHEVT( path, strlen(path), idn, 0, 2);
#endif

		if ( PAWCHN.nchevt == 0 ) {
			char	buf[81];
			PFSOCK( 0 );
			sprintf( buf, "HGETNT:%s", id_string );
			buf[80] = '\0';

			CZPUTA( buf, ierr );
			HLOGPF( " ", ierr );
#ifdef LONG_LONG_CHAIN
			{
		        	long long int size;
		        	int one = 1;
		        	int len;
		        	len = strlen(path);
		        	size = QUEST.iquest[0];
				PAWCHN.nchevt = pchevt_(
					path, &len, &idn, 
					&size, &one );
			}
#else
			PAWCHN.nchevt = PCHEVT(
				path, strlen(path), idn, 
				QUEST.iquest[0], 1 );
#endif
		}
	} 
#if DEBUGGING_PCNEXT_OPTIMIZATION
	dump_chain_table(log);
	fclose(log);
#endif

        if ( ! PNTOLD.ntold )  HBNAME( idn, " ", izero, "$CLEAR" );

	/* done */

	return R_NOERR;
}

#if defined PCNEXT_OPTIMIZATION

/*  The following functions replace equivalent (more or less) counterparts
 *  in CERNLIB:
 *	h_range_piaf() 	- replaces subroutine HRNGPF() in packlib/hbook/code
 *	h_next_nt()	- replaces subroutine PCNEXT() in pawlib/paw/code
 *
 *  The standard library routines implement a simple partitioning scheme
 *  for dividing the work of ntuple processing between the slaves on piaf.
 *  Under the standard scheme, a ntuple or ntuple chain is split into N
 *  equal segments and the segment S is assigned to slave S, for S=1..N.
 *
 *  One improvement to this scheme has already been implemented, under the
 *  compilation option WEIGHTED_NTUPLE_PARTITIONING, which allows the
 *  user to specify relative weight factors to the slaves to direct the
 *  relative processing load assigned to each.  The weights are loaded
 *  using a custom message from the piaf client, "piaf/message reload s w".
 *
 *  The basic piaf design requires that the ntuple data be served to all
 *  nodes from a common network file system, so that all slaves can
 *  simultaneously access all elements of a chain.  This is used, for
 *  example, to allow all slaves to initialize an ntuple chain by opening
 *  all of the files in the chain at the time the ntuple is loaded.  But
 *  once the ntuple is loaded, the slaves only read data from the files
 *  that are a part of their assigned segment, and performance is greatly
 *  improved if it can be arranged that the data required by each slave
 *  are actually found in files being served from its local file system.
 *  One way to accomplish this is to pre-distribute the files across
 *  the piaf slaves such that the simple piaf partitioning scheme makes
 *  a slave's segment overlap well with its local file set.  However, it
 *  makes much more sense to design a more flexible piaf partitioning
 *  algorithm that does not require a rigid arrangement of data across
 *  multiple nodes in a network.
 *
 *  The new piaf partitioning scheme is enabled with the obscure compiler
 *  option PCNEXT_OPTIMIZATION.  To take advantage of this feature, the
 *  user must specify the filenames in a ntuple chain beginning with the
 *  prefix //piaf and ending with the suffix :S where S is the index of
 *  the piaf slave that is supposed to access the data stored in that
 *  file.  Any files in the chain which do not carry the :S suffix will
 *  be split up among the slaves according to the old agorithm, so the
 *  new scheme is fully backward-compatible with the standard one, with
 *  or without WEIGHTED_NTUPLE_PARTITIONING.  Mixing of suffixed and plain
 *  chain members is permitted.  The following rules apply.
 *
 *   1. The standard bounds (see paw command nt/plot) "nevent", "ifirst"
 * 	are still respected.  Events outside the specified bounds are
 *	ignored, whether or not they fall inside suffixed chain files.
 *
 *   2. Weighting of slave processing loads, if selected, only applies to
 *	the processing of unsuffixed chain files.  If a suffix is specified
 *	then it overrides any claim that another slave may think it has to
 *	process the data in a file.
 *
 *   3. In contrast to the behaviour of HRNGPF, h_range_piaf may return
 *	a limited segment of the complete range of ntuple rows that the slave
 *	is supposed to process.  It does this because under the new scheme
 *	it frequently happens that the events a given slave is supposed to
 *	process do not appear in a single contiguous segment.  In that case,
 *	it is up h_next_nt to update the lower and upper bounds on the range
 *	to be processed each time a new file in a chain is loaded.  It is up
 *	to the query processor to notice those updates and respect them.
 *
 *   4. The new scheme obviously has benefits only in the case of a ntuple
 *	chain, where more than one ntuple file is being processed.  The
 *	:S suffix may also be used on individual-file ntuples on piaf, but
 *	it is ignored and the old-style partitioning is used.
 *
 *   5. If the user makes a mistake and specifies suffix :S with a value of
 *	S that is less than 1 or greater than the number of slaves then it
 *	is ignored and treated as if the file were unsuffixed.
 */

int
h_range_piaf(
	int		slave_count,
	int		slave_index,
# ifdef LONG_LONG_CHAIN
	long long int *	first_event,
	long long int *	final_event
# else
	int *		first_event,
	int *		final_event
# endif
) {
	int	file;
	int	slice;
	int	slice_count = slave_count +1;
	long long int	slice_rows[slice_count];
	long long int	shared_rows[slice_count];
	long long int	skip_count, use_count;
	long long int	row, total_rows;

# if DEBUGGING_PCNEXT_OPTIMIZATION
	int d;
	FILE *log;
	char filename[99];
	sprintf(filename,"piafs%d.dlog",slave_index);
	log=fopen(filename,"a");
	fprintf(log,"Entry to h_range_piaf( ");
	fprintf(log,"slave_count = %d, ",slave_count);
	fprintf(log,"slave_index = %d, ",slave_index);
#  ifdef LONG_LONG_CHAIN
	fprintf(log,"first_event = %lld, ",*first_event);
	fprintf(log,"final_event = %lld )\n",*final_event);
#  else
	fprintf(log,"first_event = %d, ",*first_event);
	fprintf(log,"final_event = %d )\n",*final_event);
#  endif
# endif
	if (chain_table_len) {
	/* compute the initial weights for each slice */
	    int slice_weight[slice_count];
	    long long int resummed_weight;
	    long long int summed_weight;
	    summed_weight = 0;
	    slice_weight[0] = 0;
            for (slice = 1; slice < slice_count; slice++) {
# ifdef WEIGHTED_NTUPLE_PARTITIONING
		int getntweight_(int *s);
		slice_weight[slice] = getntweight_(&slice);
# else
		slice_weight[slice] = 1;
# endif
		summed_weight += slice_weight[slice];
	    }
# if DEBUGGING_PCNEXT_OPTIMIZATION
	    fprintf(log,"Initial chain table :\n");
	    for (d = 0; d < chain_table_len; d++) {
		fprintf(log," %4d : %10d %2d %10d %10d\n",d,
			chain_table[d].nrows,
			chain_table[d].slice,
			chain_table[d].start,
			chain_table[d].end);
	    }
# endif
	/* count the number of rows preallocated to each slice */
	    for (slice = 0; slice < slice_count; slice++) {
		slice_rows[slice] = 0;
	    }
	    total_rows = row = 0;
	    for (file = 0; file < chain_table_len; file++) {
		int nrows = chain_table[file].nrows;
		int from, to;
		if (*first_event > row+nrows || *final_event <= row) {
		    chain_table[file].start = 1;
		    chain_table[file].end = 0;
		    row += nrows;
		    continue;
		}
		from = max(*first_event-row,1);
		to = min(*final_event-row,nrows);
		slice = chain_table[file].slice;
		slice_rows[slice] += to-from+1;
		total_rows += to-from+1;
		chain_table[file].start = from;
		chain_table[file].end = to;
		row += nrows;
	    }
# if DEBUGGING_PCNEXT_OPTIMIZATION
	    fprintf(log,"Initial work table :\n");
	    for (d = 0; d < slice_count; d++) {
		fprintf(log," %4d : %10d %15lld\n",d,
			slice_weight[d],
			slice_rows[d]);
	    }
	    fprintf(log," total: %10lld %15lld\n",
			summed_weight,
			total_rows);
# endif
	/* recompute the weights based on preallocated work */
	    resummed_weight = 0;
            for (slice = 1; slice < slice_count; slice++) {
		if (slice_rows[slice]*summed_weight > 
			total_rows*slice_weight[slice]) {
                 /*
                  * This slave is too busy to be allowed to
                  * contribute to the work in slice 0, but he
                  * must still complete his own slice_rows[slice].
                  */
		    slice_weight[slice] = 0;
		}
		resummed_weight += slice_weight[slice];
	    }
# if DEBUGGING_PCNEXT_OPTIMIZATION
	    fprintf(log,"Final work table :\n");
	    for (d = 0; d < slice_count; d++) {
		fprintf(log," %4d : %10d %15lld\n",d,
			slice_weight[d],
			slice_rows[d]);
	    }
	    fprintf(log," total: %10lld %15lld\n",
			resummed_weight,
			total_rows);
# endif
	/* allocate a fair share of slice 0 for this slave */
	    summed_weight = 0;
	    shared_rows[0] = 0;
            for (slice = 1; slice < slice_count; slice++) {
		summed_weight += slice_weight[slice];
		shared_rows[slice] = slice_rows[0]*summed_weight
						/resummed_weight;
	    }
	    skip_count = shared_rows[slave_index-1];
	    use_count = shared_rows[slave_index] - skip_count;
	    for (file = 0; file < chain_table_len; file++) {
		int from = chain_table[file].start;
		int to = chain_table[file].end;
		if (chain_table[file].slice > 0) {
		    if (chain_table[file].slice != slave_index) {
			chain_table[file].end = 0;
		    }
		    continue;
		}
		int nrows = to-from+1;
		if (skip_count >= nrows) {
		    chain_table[file].end = 0;
		    skip_count -= nrows;
		    continue;
		}
		chain_table[file].start += skip_count;
		nrows -= skip_count;
		skip_count = 0;
		use_count -= nrows;
		if (use_count < 0) {
		    chain_table[file].end += use_count;
		    use_count = 0;
		}
	    }
	/* find the first non-empty segment for this slave and return it */
	    row = 0;
	    *final_event = 0;
	    for (file = 0; file < chain_table_len; file++) {
		int from = chain_table[file].start;
		int to = chain_table[file].end;
		if (to >= from) {
		    *first_event = row+from;
		    *final_event = row+to;
		    break;
		}
		row += chain_table[file].nrows;
	    }
# if DEBUGGING_PCNEXT_OPTIMIZATION
	    fprintf(log,"Final chain table :\n");
	    for (d = 0; d < chain_table_len; d++) {
		fprintf(log," %4d : %10d %2d %10d %10d\n",d,
			chain_table[d].nrows,
			chain_table[d].slice,
			chain_table[d].start,
			chain_table[d].end);
	    }
# endif
	}

	/* otherwise, just apply the standard partitioning algorithm */

	else {
# ifdef WEIGHTED_NTUPLE_PARTITIONING
	    long long int preWeight=0;
	    long long int postWeight=0;
	    long long int totalWeight=0;
	    int getntweight_(int *s);
	    int is;
	    for (is=1; is <= slave_count; is++) {
		if (is == slave_index) {
		    preWeight = totalWeight;
		    totalWeight += getntweight_(&is);
		    postWeight = totalWeight;
		}
		else {
		    totalWeight += getntweight_(&is);
		}
	    }
	    total_rows = *final_event - *first_event +1;
	    *final_event = *first_event + total_rows*postWeight/totalWeight -1;
	    *first_event = *first_event + total_rows*preWeight/totalWeight;
# else
	    total_rows = *final_event - *first_event +1;
	    *final_event = *first_event + total_rows*slave_index/slave_count -1;
	    *first_event = *first_event + total_rows*(slave_index-1)/slave_count;
# endif
	}
#if DEBUGGING_PCNEXT_OPTIMIZATION
	fprintf(log,"Return from h_range_piaf( ");
	fprintf(log,"slave_count = %d, ",slave_count);
	fprintf(log,"slave_index = %d, ",slave_index);
# if LONG_LONG_CHAIN
	fprintf(log,"first_event = %lld, ",*first_event);
	fprintf(log,"final_event = %lld );\n",*final_event);
# else
	fprintf(log,"first_event = %d, ",*first_event);
	fprintf(log,"final_event = %d );\n",*final_event);
# endif
	fclose(log);
#endif
}

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
	int *	 	next_evt,
	int *	 	end_evt
# endif
) {
	int	nextfile;
# ifdef LONG_LONG_CHAIN
	long long int nextrow = 1;
	long long int ndone = *nchrow + *nrow;
# else
	int 	nextrow = 1;
	int	ndone = *nchrow + *nrow;
# endif
# if DEBUGGING_PCNEXT_OPTIMIZATION
	FILE *log;
	char filename[99];
	sprintf(filename,"piafs%d.dlog",HCPIAF.mysid);
	log=fopen(filename,"a");
	fprintf(log,"Entry to h_next_nt( ");
	fprintf(log,"id = %d, ",id);
	fprintf(log,"nchrow = %lld, ",*nchrow);
	fprintf(log,"ndim = %d, ",*ndim);
	fprintf(log,"nrow = %d, ",*nrow);
	fprintf(log,"iend = %d, ",*iend);
	fprintf(log,"next_evt = %lld, ",*next_evt);
	fprintf(log,"end_evt = %lld )\n",*end_evt);
# endif
	for (nextfile = 0; nextfile < chain_table_len;) {
		nextrow += chain_table[nextfile++].nrows;
		if (nextrow > ndone)
			break;
	}
	while (nextfile < chain_table_len && 
	       chain_table[nextfile].start > chain_table[nextfile].end) {
		char chain[80];
		char tentry[80];
		int lc,lt;
		PCHNXT(chain,lc,tentry,lt);
		*nchrow += chain_table[nextfile++].nrows;
	}
# ifdef LONG_LONG_CHAIN
	pcnext_(&id,nchrow,ndim,nrow,iend );
# else
	PCNEXT(id,*nchrow,*ndim,*nrow,*iend);
# endif
	*iend = (nextfile == chain_table_len);
	*nrow = chain_table[nextfile].nrows;
        if (chain_table[nextfile].start <= chain_table[nextfile].end) {
	    *next_evt = *nchrow + chain_table[nextfile].start;
	    *end_evt = *nchrow + chain_table[nextfile].end;
	}
# if DEBUGGING_PCNEXT_OPTIMIZATION
	fprintf(log,"Return from h_next_nt( ");
	fprintf(log,"id = %d, ",id);
	fprintf(log,"nchrow = %lld, ",*nchrow);
	fprintf(log,"ndim = %d, ",*ndim);
	fprintf(log,"nrow = %d, ",*nrow);
	fprintf(log,"iend = %d, ",*iend);
	fprintf(log,"next_evt = %lld, ",*next_evt);
	fprintf(log,"end_evt = %lld )\n",*end_evt);
	fclose(log);
# endif
}

int h_nrows_piaf ()
{
	int nrows, c;
	for (c=0; c < chain_table_len; c++) {
		nrows += chain_table[c].end-chain_table[c].start+1;
		if (nrows < 0) {
			return -1;
		}
	}
	return nrows;
}
#endif

static void
hntld_C( char * name )
{
	int	dummy_id;
	char	*dummy_path;

	dummy_path = 0;

	h_load_nt( name, &dummy_path, &dummy_id );

	if (dummy_path) str_del( dummy_path );

	return;
}

FCALLSCSUB1(hntld_C,HNTLD,hntld,STRING)


static void
hrstd_C( )
{
/* Restore the old path stored by h_load_nt from fortran */

	HCDIR( old_path, " " );

	return;
}

FCALLSCSUB0(hrstd_C,HRSTD,hrstd)

void
h_reset_dir( )
{
/* Restore the old path stored by h_load_nt from C */

	HCDIR( old_path, " " );

	return;
}


int
h_load_histo( char *id_string, int * idp, int *id_dim )
{
	char	*p;
	long	id;
	int	err = R_NOERR;

	id = strtol( id_string, &p, 10 );

	if ( (strlen(id_string) == 0) || (*p != 0) ) {
		return R_SYNTAX_ERROR;
	}

	HFIND( id, "h_load_histo" );
	if ( HCBOOK.lcid == 0 ) {
		return R_NAME_UNDEF;
	}

	HDCOFL();

	if ( (HCBITS.i1 != 1) && (HCBITS.i2 !=1) && (HCBITS.i8 != 1) ) {
		return R_NO_HISTOGRAM;
	}

	*idp = id;

	if ( (HCBITS.i1 == 1) || (HCBITS.i8 == 1) ) {
		*id_dim = 1;
	} else {
		*id_dim = 2;
	}

	return err;
}


/*
 *  Return the number of rows in a ntuple or chain
 */
#ifdef LONG_LONG_CHAIN
long long int
#else
void
#endif
h_hnoent( int idn, bool use_chain, int * ep )
{
	int	entries;

	if ( use_chain && PAWCHN.chain ) {
		*ep = PAWCHN.nchevt;
#ifdef LONG_LONG_CHAIN
		return PAWCHN.nchevt;
#endif
	} else {
		HNOENT( idn, entries );
		*ep = entries;
#ifdef LONG_LONG_CHAIN
		return entries;
#endif
	}
}


/*
 *  Return the number of columns in a ntuple
 */
void
h_hnocol( int idn, int * colp )
{
	HFIND(idn,"HNOENT");

	qp_assert( HCBOOK.lcid != 0); /* found ntuple */
	qp_assert((PAWC.iq_eqv[HCBOOK.lcid+KBITS-1]&0x8)==0x8); /* is ntuple */

	*colp = PAWC.iq_eqv[HCBOOK.lcid+ZNDIM-1];
}


bool
h_flag_1d(
	int		id
){
	HFIND(id,"h_flag_1d");

	if ( HCBOOK.lcid == 0 ) {
		sf_report( "h_flag_1d: Id %d not found\n", id );
		return FALSE;		/* id not found */
	}

	if ( (PAWC.iq_eqv[HCBOOK.lcid+KBITS-1] & 0x1) == 0x1 ) {
		return TRUE;
	} else {
		return FALSE;
	}
}


bool
h_flag_2d(
	int		id
){
	HFIND(id,"h_flag_2d");

	if ( HCBOOK.lcid == 0 ) {
		sf_report( "h_flag_2d: Id %d not found\n", id );
		return FALSE;		/* id not found */
	}

	if ( (PAWC.iq_eqv[HCBOOK.lcid+KBITS-1] & 0x2) == 0x2 ) {
		return TRUE;
	} else {
		return FALSE;
	}
}


bool
h_flag_profile(
	int		id
){
	HFIND(id,"h_flag_profile");

	if ( HCBOOK.lcid == 0 ) {
		sf_report( "h_flag_profile: Id %d not found\n", id );
		return FALSE;		/* id not found */
	}

	if ( (PAWC.iq_eqv[HCBOOK.lcid+KBITS-1] & 0x80) == 0x80 ) {
		return TRUE;
	} else {
		return FALSE;
	}
}


char	h_rwn_chtitl[128];
int	h_rwn_nvar;
char	h_rwn_tags[MAX_RWN_COLS][9];
float	h_rwn_rlow[MAX_RWN_COLS];
float	h_rwn_rhigh[MAX_RWN_COLS];


bool
h_rwn_getInfo(
	int		idn
) {
	int		i;
	char		*p;

	h_rwn_nvar = MAX_RWN_COLS;
	HGIVEN(idn,h_rwn_chtitl,h_rwn_nvar,h_rwn_tags,h_rwn_rlow,h_rwn_rhigh);

	if ( h_rwn_nvar == 0 ) {
		return FALSE;
	}

	for( i=0 ; i < h_rwn_nvar ; i++ ) {
		for ( p = h_rwn_tags[i]; *p == ' ' ; p++ ) { }
		memmove( h_rwn_tags[i], p, 9 - ( p - h_rwn_tags[i] ) );
	}

	return TRUE;
}


int
h_rwn_getIndex(
	int		idn,
	char * const	name
) {
	int		nvar, i;

	if ( ! h_rwn_getInfo( idn ) ) {
		return -1;
	}

	for( i=0 ; i < h_rwn_nvar ; i++ ) {

		if ( strcasecmp( name, h_rwn_tags[i] ) == 0 ) {
			return i;
		}
	}

	return -1;
}


SMap
h_get_labels(
	int		id,
	char		*chopt
) {
	SMap	sm;
	int	i, j, n;
	char	lab[MAX_LABEL_LEN+1];

	qp_assert( !strcmp(chopt,"X") ||
		!strcmp(chopt,"Y") || !strcmp(chopt,"Z") );

	n = HLABNB(id,chopt);
	sm = smap_new( n );

	for ( i=0 ; i < n ; i++ ) {
		lab[0] = '\0';
		HLGNXT(id,i+1,lab,chopt);
		j = strlen( lab );
		if ( j < MAX_LABEL_LEN ) {
			memset( &lab[j], ' ', MAX_LABEL_LEN - j );
			lab[MAX_LABEL_LEN] = '\0';
		}
		smap_add( sm, lab );
	}

	return sm;
}


void
h_hbook1_labels(
	int		idh,
	char *		title,
	SMap		label_list
)
{
	int	n, i;
	char	clab[MAX_MAP_STRINGS][MAX_LABEL_LEN+1];

	for( i=0 ; i < MAX_MAP_STRINGS ; i++ ) {
		clab[i][0] = '\0';
	}

	n = smap_entries(label_list);
	HBOOK1( idh, title, n, 1., (float) (n+1), 0. );

	for ( i=0 ; i < n ; i++ ) {
		strncpy( clab[i], smap_get( label_list, i ), MAX_LABEL_LEN);
		clab[i][MAX_LABEL_LEN] = '\0';
	}
	HLABEL( idh, n, clab, "NX" );
}


void
h_hbook2_labels(
	int		idh,
	char *		title,
	SMap		labelx_list,
	SMap		labely_list,
	int		nbin[],
	float		min[],
	float		max[]
)
{
	int	n, i;
	char	clab[MAX_MAP_STRINGS][MAX_LABEL_LEN+1];

	for( i=0 ; i < MAX_MAP_STRINGS ; i++ ) {
		clab[i][0] = '\0';
	}

	if ( labely_list != 0 ) {
		nbin[0] = smap_entries(labely_list);
		min[0] = 1.;
		max[0] = nbin[0] + 1.;
	}

	if ( labelx_list != 0 ) {
		nbin[1] = smap_entries(labelx_list);
		min[1] = 1.;
		max[1] = nbin[1] + 1.;
	}

	HBOOK2( idh, title,
		nbin[1], min[1], max[1],
		nbin[0], min[0], max[0],
		0. );

	if ( labely_list != 0 ) {
		for ( i=0 ; i < nbin[0] ; i++ ) {
			strncpy( clab[i], smap_get( labely_list, i ),
				MAX_LABEL_LEN);
			clab[i][MAX_LABEL_LEN] = '\0';
		}
		HLABEL( idh, nbin[0], clab, "NY" );
	}

	if ( labelx_list != 0 ) {
		for ( i=0 ; i < nbin[1] ; i++ ) {
			strncpy( clab[i], smap_get( labelx_list, i ),
				MAX_LABEL_LEN);
			clab[i][MAX_LABEL_LEN] = '\0';
		}
		HLABEL( idh, nbin[1], clab, "NX" );
	}
}

#if DEBUGGING_PCNEXT_OPTIMIZATION
int dump_chain_table (FILE* log)
{
        int i;
	int nrows;
	fprintf(log,"dump of chain_table:\n");
	for (nrows=i=0; i < chain_table_len; i++)
	{
		fprintf(log,"  nrows=%d, slice=%d, start=%d, end=%d\n",
			chain_table[i].nrows,
			chain_table[i].slice,
			chain_table[i].start,
			chain_table[i].end);
		nrows += chain_table[i].nrows;
        }
        fprintf(log,"total rows %d\n\n",nrows);
	return 0;
}
#endif
