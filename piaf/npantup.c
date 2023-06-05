/* 
 * -R.T. Jones  April 23, 2005
 *  Optimization of PCNEXT processing introduced for parallel
 *  processing on a large number of slaves.  It avoids some
 *  significant overhead associated with open/close cycles on
 *  files that are never accessed.   For more information, see
 *  qp_hbook_if.c source code.  To disable this feature, comment
 *  out the following define below and in qp_hbook_if.c source.
 */
#define PCNEXT_OPTIMIZATION 1
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
 * -R.T. Jones  July 15, 2002
 *  Added support for the progress meter for piaf ntuple chain
 *  processing.  It used to work with the original piaf service
 *  but seems to never have been implemented for the new qp.
 *  To disable this feature, comment out the following line.
 */
#define PIAF_PROGRESS_METER 1
/*
 * -R.T. Jones, July 12, 2002
 *  Added a couple of fixes to make things work better under Linux.
 *  To disable these changes, comment out the following line.
 */
#define MAKE_GCC_HAPPY 1
/*
 *  npantup.c  --
 *	Kuip action routine for ntuple/{loop,project,plot}
 *
 *  Original: 26-Mar-1995 18:37
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: npantup.c,v 1.4 2006/03/30 02:22:04 radphi Exp $
 *
 *  $Log: npantup.c,v $
 *  Revision 1.4  2006/03/30 02:22:04  radphi
 *  npantup.c, qp_execute.c
 *     - changed back to report statistics without rescaling -rtj
 *  qp_hbook_if.c
 *     - bug fix to handle case with only one event per file -rtj
 *     - bug fix to handle switching between ntuples within a given chain -rtj
 *     - other minor fixes -rtj
 *  pfpaw.F
 *     - fixed handling of work dialog boxes to prevent repeated popping -rtj
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
 *  Revision 1.2  2004/07/01 17:59:20  radphi
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
 *  Revision 1.95  2001/09/25 12:44:59  couet
 *  - Option S in the NT/PLOT command, used the limit of the histogram created
 *    by the previous NT/PLOT command. Example:
 *
 *    nt/plot 10.age age<40
 *    null 0 100 0 100
 *    nt/plot 10.age
 *
 *    The last NT/PLOT command was also done with the cut age<40
 *
 *  Revision 1.94  2001/09/18 13:41:35  couet
 *  - In case of alphanumeric label on 1D or 2D histograms, the alphabetical
 *    sorting done by smap_sort was wrong. With the ntuple generated by the
 *    following example ...
 *
 *        COMMON /PAWC/  SPACE(50000)
 *
 *        CHARACTER*40   CHFILE
 *        CHARACTER*4    TYPE
 *        COMMON /CEXAM/ TYPE
 *
 *        CALL HLIMIT ( 50000 )
 *
 *        CHFILE = 'test.cwn'
 *        CALL HROPEN ( 11, 'MYFILE', CHFILE, 'N', 1024, ISTAT )
 *        CALL HBNT   ( 10, 'TEST NTUPLE', ' ' )
 *        CALL HBNAMC ( 10, 'TYPE', TYPE, 'TYPE:C' )
 *
 *        CALL FILL('BB',1)
 *        CALL FILL('CC',2)
 *        CALL FILL('EE',3)
 *        CALL FILL('FF',4)
 *        CALL FILL('GG',5)
 *        CALL FILL('HH',6)
 *        CALL FILL('IT',7)
 *        CALL FILL('JJ',8)
 *        CALL FILL('LL',9)
 *        CALL FILL('MM',8)
 *        CALL FILL('NN',7)
 *        CALL FILL('OO',6)
 *        CALL FILL('SE',5)
 *        CALL FILL('ZZ',4)
 *        CALL FILL('AA',3)
 *
 *        CALL HROUT  ( 10, ICYCLE, ' ' )
 *        CALL HREND  ( 'MYFILE' )
 *        END
 *
 *        SUBROUTINE FILL (C,N)
 *        CHARACTER*(*) C
 *        CHARACTER*4    TYPE
 *        COMMON /CEXAM/ TYPE
 *        TYPE = C
 *        DO I=1,N
 *           CALL HFNT(10)
 *        ENDDO
 *        END
 *
 *    ... the entries in the SE channel were dropped into the underflow channel
 *    when one use the command:
 *
 *    nt/plot 10.type
 *
 *    The smap_sort sorting has been replaced by the sorting done by the HBOOK
 *    routine HLABEL.
 *
 *  Revision 1.93  1999/08/31 08:47:33  couet
 *  - A new set of PAW commands using Neural Networks.
 *    These new commands are using the MLPfit package.
 *
 *  Revision 1.92  1999/07/05 15:43:31  couet
 *  - hbook_interface.h in now replaced by hbook.h in CVSCOSRC
 *
 *  Revision 1.91  1999/07/02 12:31:45  couet
 *  - hplot_interface.h is now replaced by hplot.h in CVSCOSRC
 *  - bug fixed in qp_execute.c: the zones were not properly reset after a
 *    NT/SCAN with option SPIDER
 *
 *  Revision 1.90  1999/07/01 14:46:59  couet
 *  - cspack_interface.h is replaced by cspack.h in CVSCOSRC
 *
 *  Revision 1.89  1999/07/01 14:18:42  couet
 *  - higz_interface.h in now replaced by higz.h in CVSCOSRC
 *
 *  Revision 1.88  1999/06/28 15:08:51  couet
 *  - use now cfortran.h in $CVSCOSRC
 *
 *  Revision 1.87  1999/03/25 14:22:07  couet
 *  - The following commands may produce wrong plot with GCUT:
 *
 *    nt/plot 1.y%z $1
 *    gcut 2 1.y%x $1 option=S
 *
 *    The S option forced GCUT to produce a histogram scatter plot
 *    which completly obscured the picture.
 *
 *  Revision 1.86  1997/02/12 14:10:43  couet
 *  if ((!cmd->u.minmax.has_events)&&(*errp == R_NOERR)){
 *
 *   should be
 *
 *  if ((!cmd_tmp->u.minmax.has_events)&&(*errp == R_NOERR)){
 *
 *   in do_gouraud
 *
 *  Revision 1.85  1996/12/12 15:21:37  lecointe
 *  Corrected a bug in expand_var_list du to a misinterpreted return value
 *  of ku_getl
 *
 *  Revision 1.84  1996/12/12 14:30:33  dinofm
 *  Code has been modified fro PIAF. Whenever an error is detected in
 *  the MINMAX command execution (e.g. Array Bound Error) a PIAF slave
 *  sends a 'Empty Histogram' synchronization sequence. This avoids
 *  all other slaves to be stuck if one of them detects such an error.
 *  All PIAF's 'Empty Histogram' synchronization sequence have been
 *  grouped in a set of functions.
 *
 *  Revision 1.83  1996/12/05 10:04:21  lecointe
 *  Two completely different things
 *  	. add a _C to "paqcut" in "qp_cuts.c" to avoid duplicate symbol
 *  on VMS
 *  	. Modified "expand_var_list" called by "do_nt_scan" in npantup.c.
 *  Removed the hand-made parser and use "ku_getl" to parse the list of
 *  variables to scan.
 *
 *  Revision 1.82  1996/10/18 13:54:49  couet
 *  - more prototype
 *
 *  Revision 1.80  1996/10/16 10:11:24  couet
 *  - mistake in a mod in do_prof_1d
 *
 *  Revision 1.79  1996/09/24 08:58:07  lecointe
 *  Add option 'N' in NTUPLE/PLOT : Fill the 1D or 2D histogram without drawing
 *  anything
 *
 *  Revision 1.78  1996/09/19 09:24:38  couet
 *  - The statistics were missing in case of 2D scatter plots (opt stat)
 *
 *  Revision 1.77  1996/09/18 14:29:33  dinofm
 *  Slaves synchronize themselves in function do_scat_2d so that the PAPLOT
 *  command is executed before any attempt to send the point to the client
 *  by IPMID. Synchronization is allowed by the PSYNC message.
 *
 *  Revision 1.76  1996/09/12 09:28:53  couet
 *  - after NTUPLE commands using a full path name ie:
 *
 *    NT/PLOT //lun1/10.x
 *
 *    The current working directory was changed to //lun1
 *
 *  Revision 1.75  1996/08/30 14:27:08  lecointe
 *  Restored the S option for 3d and 4d scatter plot
 *
 *  Revision 1.74  1996/08/30 10:04:25  lecointe
 *  Restored Gouraud Shading in NT/PLOT
 *
 *  Revision 1.73  1996/08/21 12:55:31  lecointe
 *  Restore the spider plot in ntuple/scan
 *
 *  Revision 1.72  1996/07/11 15:00:08  couet
 *  - The error messages printed in calc_event_range in case of a wrong number
 *    of events are improved. In particular for empty Ntuples or ifirst greater
 *    than the number of events, the message was wrong.
 *
 *  Revision 1.71  1996/07/04 13:25:55  couet
 *  - The LOCATOR on ntuples (in PAW and PAW++) was not working.
 *    The calls to IGPID were missing.
 *
 *  Revision 1.70  1996/06/05 08:09:01  couet
 *  - The user title (UPT UTIT and TITLE with option U) was not taken
 *    into account in the NT/PLOT command
 *
 *  Revision 1.69  1996/05/24 09:14:21  dinofm
 *  Memory leakage on id_path fixed
 *
 *  Revision 1.68  1996/05/22 09:44:54  dinofm
 *  When an histogram id is specified in NT/PL, the related histogram is
 *  cleaned by HRESET before being pushed.
 *
 *  Revision 1.67  1996/05/21 15:59:43  couet
 *  - It was not allowed to do the "S"ame option in NT/PLOT if the previous
 *    plot was not done with NT/PLOT. Now, in that case, the coordinates for
 *    the Same are taken from IGQWK.
 *
 *  Revision 1.66  1996/05/15 13:11:36  maartenb
 *  - Fix the CSELECT command.
 *
 *  Revision 1.65  1996/05/13 16:22:06  dinofm
 *  NT/LOOP of Piaf ntuples is now allowed.
 *
 *  Revision 1.64  1996/05/10 12:23:30  dinofm
 *  PIAF debug messages (Ntuple on Piaf etc.) are issued only if the 'verbose'
 *  flag is not 0.
 *
 *  Revision 1.63  1996/05/09 10:17:00  dinofm
 *  Fixed 'No event selected' code for NT/PL & NT/GCUT. Fixed string allocation
 *  when 'old selection function' detected on PIAF.
 *
 *  Revision 1.62  1996/05/08 16:13:03  dinofm
 *  NT/SCAN on PIAF is enabled again.
 *
 *  Revision 1.61  1996/05/06 13:34:16  dinofm
 *  Code modified to take care of empty histograms detection on slave(s).
 *
 *  Revision 1.60  1996/04/30 10:10:50  maartenb
 *  - Fix iand cleanup of the detection and handling of profile histograms.
 *
 *  Revision 1.59  1996/04/24 08:24:02  dinofm
 *  The warning message for obsolete NT/PLOT default selection function
 *  has been removed.
 *
 *  Revision 1.58  1996/04/23 18:38:07  maartenb
 *  - Add RCS keywords
 *
 *
 */

#include	<string.h>
#include	<ctype.h>
#include	<float.h>
#include	<math.h>
#include	<errno.h>

#include	"str.h"		/* IRIX stdarg.h clash ... */

#ifdef MAKE_GCC_HAPPY
# ifdef linux
# define f2cFortran 1
# endif
#endif
#include	<cfortran/cfortran.h>

#define	CF_TRUE		C2FLOGICAL(1)
#define	CF_FALSE	C2FLOGICAL(0)

#include	"c_hcpiaf.h"
#include	"c_pawchn.h"
#include	"c_pcaddr.h"
#include	"c_pccsel.h"
#include	"c_pntold.h"
#include	"c_quest.h"
#ifdef PIAF_PROGRESS_METER
#include	"c_pcntwk.h"
#endif


#include	"errors.h"
#include	"hash_int_table.h"
#include	"hbook_defs.h"
#include	<cfortran/hbook.h>
#include	<cfortran/higz.h>
#include	<cfortran/hplot.h>
#include	"kuip_interface.h"
#include	"paw_interface.h"
#include	<cfortran/cspack.h>
#include	"qp_command.h"
#include	"qp_compile.h"
#include	"qp_cuts.h"
#include	"qp_execute.h"
#include	"qp_hbook_if.h"
#include	"qp_plot_opt.h"
#include	"qp_query.h"
#include	"qp_report.h"
#include	"qpflags.h"
#include	"svec.h"

extern void PiafNoEvt ( bool );
extern void PackLabels( SMap, int *, String, int * );
extern void UnpackLabels( SMap, String, int );
extern int NewPiaf();

float text_angle = 0.;  /* Angle used for comments printed with IGTEXT */
/*
 *  The following structures are used to make option S working
 */

typedef struct _zone_hist_ {
	int		idim;		/* dimensionality of the zone */
	int		nxbin;		/* # of x bins */
	float		xlow, xhigh;	/* booking limits in x */
	int		nybin;		/* # of y bins */
	float		ylow, yhigh;	/* booking limits in y */
} ZoneHist;

static ZoneHist	theLastZone;


static void
PiafEmpty1D (char *title, bool labels) {
  PFEMPTY ( title );
  PiafNoEvt( labels );
}
  
static void
PiafEmpty2D (char *title, bool labelx, bool labely) {
  int i;
  PFEMPTY( title );
  for( i=0 ; i < 2 ; i++ ) {
    if ( i==0 && labely ) {
      PiafNoEvt( TRUE );
      continue;
    }
    if ( i==1 && labelx ) {
      PiafNoEvt( TRUE );
      continue;
    }
    PiafNoEvt( FALSE );
  }
}

static void
PiafEmptyNoLabels (char *title, int dimension) {
  int i;
  PFEMPTY( title );
  for( i=0 ; i < dimension ; i++ ) {
    PiafNoEvt( FALSE );
  }
}

static void
#ifdef LONG_LONG_CHAIN
calc_event_range( int id, long long int *firstp, long long int *neventp )
{
	long long int	entries, ifirst, nevent, tmp;
#else
calc_event_range( int id, int *firstp, int *neventp )
{
	int	entries, ifirst, nevent, tmp;
#endif
	int ent;

	ifirst = *firstp;
	nevent = *neventp;

	entries = h_hnoent( id, TRUE, &ent );

	if ( entries == 0 ) {
		sf_report( "Ntuple is empty.\n");
		nevent = 0;
		ifirst = 1;
	} else if ( ifirst > entries ) {
#ifdef LONG_LONG_CHAIN
		sf_report( "Ntuple has only %lld events. "
#else
		sf_report( "Ntuple has only %d events. "
#endif
			"No events processed.\n",
			entries);
		nevent = 0;
		ifirst = 1;
	} else if ( nevent == 99999999 ) {
	        nevent = entries - ifirst + 1;
	} else if ( ifirst + nevent - 1 > entries ) {
		tmp = entries - ifirst + 1;
#ifdef LONG_LONG_CHAIN
		sf_report( "Ntuple has only %lld events. "
				"%lld events processed.\n",
#else
		sf_report( "Ntuple has only %d events. "
				"%d events processed.\n",
#endif
				entries, tmp );
		nevent = tmp;
	}

#ifdef PIAF_PROGRESS_METER
	PCNTWK.iminev = 1;
        if (PCNTWK.nchain > 0) {
                PCNTWK.imaxev = (nevent < 999999999)? nevent : 999999999;
        }
        else {
                PCNTWK.imaxev = nevent;
        }
#endif

	if ( HCPIAF.slavpf ) {
#ifdef LONG_LONG_CHAIN
		long long int ito;
#else
		int ito;
#endif

#ifdef PCNEXT_OPTIMIZATION
		ito = ifirst + nevent - 1;
		h_range_piaf(HCPIAF.ngsize, HCPIAF.mysid, &ifirst, &ito);
#elif WEIGHTED_NTUPLE_PARTITIONING
		long long int preWeight=0;
		long long int postWeight=0;
		long long int totalWeight=0;
		int getntweight_(int *s);
		int is;
		for (is=1; is <= HCPIAF.ngsize; is++) {
			if (is == HCPIAF.mysid) {
				preWeight = totalWeight;
				totalWeight += getntweight_(&is);
				postWeight = totalWeight;
			}
			else {
				totalWeight += getntweight_(&is);
			}
		}
		ito = ifirst + nevent*postWeight/totalWeight -1;
		ifirst = ifirst + nevent*preWeight/totalWeight;
#elif LONG_LONG_CHAIN
		ito = ifirst + nevent*HCPIAF.mysid/HCPIAF.ngsize -1;
		ifirst = ifirst + nevent*(HCPIAF.mysid-1)/HCPIAF.ngsize;
#else
		ito = ifirst + nevent - 1;
		HRNGPF( HCPIAF.ngsize,  HCPIAF.mysid, ifirst, ito);
#endif


		if ( ito - ifirst < 0 ) {
			nevent = 0;	/* fewer events than slaves ! */
		} else {
			nevent = ito - ifirst + 1;
		}
	}

	*firstp = ifirst;
	*neventp = nevent;
}


static int
split_id_string (
	char *		expr_string,
	String *	id_string,
	SVec *		expressions
)
{
	char	*p, *q, *s;
	int	n;
	String	es;
	SVec	sv;

	es = str_new( expr_string );

	sv = svec_new( MAX_EXPRS );
	*expressions = sv;

	/* everything before the dot is the ntuple id */
	p = strchr( es, '.' );

	if ( p == 0 ) {		/* Only an ntuple id, no expressions */
		*id_string = es;
		return R_NOERR;
	} 

	*p = '\0';
	p += 1;
	*id_string = str_new( es );

	/* a possible cycle number must be moved to id_string */
	q = strrchr( p, ';' );

	if ( q != 0 ) {
		*id_string = str_merge( *id_string, str_new( q ) );
		*q = '\0';
		q += 1;
	}

	/* split the remaining part into '%' separated expressions */

	n = 0;
	s = p;
	if ( *s != '\0' ) {
		while( s != 0 ) {
			q = strchr( s, '%' );
			if ( q != 0 ) {
				*q = '\0';
				q += 1;
			}

			if ( *s == '\0' ) {
				/* syntax error */
				sf_report(
				"Error in syntax (expr [ '%' expr [...] ])\n" );
				str_del( es );
				str_del( *id_string );
				svec_del( sv );
				return R_SYNTAX_ERROR;	/* mem leak */
			}

			if ( (q != 0) && (*q == '\0') ) {
				/* syntax error */
				sf_report(
				"Error in syntax (expr [ '%' expr [...] ])\n" );
				str_del( es );
				str_del( *id_string );
				svec_del( sv );
				return R_SYNTAX_ERROR;
			}

			if ( svec_entries(sv) == MAX_EXPRS ) {
				/* syntax error */
				sf_report(
					"To many expressions (max %d)\n",
					MAX_EXPRS );
				str_del( es );
				str_del( *id_string );
				svec_del( sv );
				return R_SYNTAX_ERROR;
			} else {
				svec_add( sv, str_new( s ) );
			}
			s = q;
		}
	}

	*expressions = sv;

	str_del( es );
	return R_NOERR;
}


#if 1
static void
calc_axis(
	int		nbin_req,
	float		min,
	float		max,
	bool		is_int,
	float *		low,
	float *		up,
	int *		nbin
) {
	double		range, base, step, tmp;

	if ( max < min ) {	/* no events */
		min = max = 0;
	}

	if ( (min < -FLT_MAX/4) || (max > FLT_MAX/4) ) {
		sf_report( "Cannot calculate proper binning for (%e,%e)\n",
			min,max);
		*low = -FLT_MAX/4;
		*up = FLT_MAX/4;
		*nbin = 64;
		return;
	}

	if ( is_int ) {
		max += 1;
	}

	/* horizontl/vertical lines */
	if ( (max - min) < 4*FLT_MIN ) {
		if ( max < 0. ) {
			max = 0.;
		} else if ( max > 0. ) {
			min = 0.;
		} else {
			max = 0.5;
			min = -0.5;
		}
		max = max * 2.;
		min = min * 2.;
	}

        range =  max - min;

	/* resolve the max -> overflow problem */
	if ( ! is_int ) {
		max += 0.001 * range;
		min -= 0.001 * range;
	}

        base = pow( 10., floor( log10( range ) ) );
        tmp = range / base;

	if ( tmp >= 5. ) {		/* 5 to 10 divisions */
            step = 1;
	} else if ( tmp >= 2. ) {	/* 4 to 10 divisions */
            step = 0.5;
	} else {			/* 5 to 10 divisions */
            step = 0.2;
	}

        tmp = min /  base;
        min = floor(tmp/step) * base * step; 
        tmp = max /  base;
        max = ceil(tmp/step)  * base * step;
        step = step * base;

	*low = min;
	*up = max;

	if ( is_int  && (step <= 10.) ) {
		step = 1.;
		*nbin = floor( (max - min) / step + 0.5);
	} else {
		nbin_req /= 10; /* maximum number of divisions == 10 */
		*nbin = floor( nbin_req * (max - min) / step + 0.5);
	}
}
#else
static void
calc_axis(
	int		nbin_req,
	float		min,
	float		max,
	bool		is_int,
	float *		low,
	float *		up,
	int *		nbin
) {
	float		dx, oldmax, oldmin, max2;
	float		tmpmin, tmpmax, bwid;
	int		nbinx;

	if ( is_int ) {
		max += 1.;
	}
	
	if ( min > max ) {
		min = -1.;
		max = 1.;
	} else if ( min == max ) {
		min -= 1.;
		max += 1.;
	} else {

		dx = (max - min ) / 100.;
		if ( dx == 0 ) {
			dx = 1.;
		}
		oldmin = min;
		oldmax = max;
		min -= 10 * dx;
		max += 10 * dx;
		if ( (oldmin >= 0.) && (min < 0.) ) {
			min = 0.1 * oldmin;
		}
		if ( (oldmax <= 0 ) && (max > 0.) ) {
			max = 0.; /* should be 0.1 *oldmax */
		}
	}

	/* code is kept identical to original :-( */
	nbinx = nbin_req;
/* init to avoid cfortran problems with uninitialised outpars */
	tmpmin = 0.; tmpmax = 0.; bwid = 0.;
	HBIN(min,max,100,tmpmin,tmpmax,nbinx,bwid);
	min = tmpmin;
	max2 = min + nbinx * bwid;
	if ( max2 <= max ) {
		nbinx = nbin_req;
		max2 = min + nbin_req * bwid;
	}
	if ( is_int ) {
		min = floor( min );
		max2 = ceil( max2 );
		nbinx = max2 - min;
		if ( nbinx > 100 ) {
			nbinx = 100;
		}
	}

	*low = min;
	*up = max2;
	*nbin = nbinx;
}
#endif


static void
add_range(
	HashIntTable	tab,
	SVec		all_tags,
	char *		var1,
	char *		var2,
	SVec		tags,
	int *		errp
) {
	int	i, i1, i2, step;

	if ( var1[0] == '\0' ) {
		i1 = 0;
	} else if ( ! HashInt_find( tab, var1, &i1 ) ) {
		sf_report( "name %s is not an ntuple variable\n", var1 );
		*errp = R_TYPE_ERROR;
		return;
	}

	if ( var2[0] == '\0' ) {
		i2 = svec_entries( all_tags ) - 1;
	} else if ( ! HashInt_find( tab, var2, &i2 ) ) {
		sf_report( "name %s is not an ntuple variable\n", var2 );
		*errp = R_TYPE_ERROR;
		return;
	}

	step = i1 <= i2 ? 1 : -1 ;
	i = i1 - step;
	do {
		i += step;
		svec_add( tags, str_new( svec_get( all_tags, i ) ) );
	} while ( i != i2 );
}


static int
expand_var_list(
	int		id,
	char *		var_list,
	SVec *		vlp
) {
	HashIntTable	tab;
	SVec		all_tags;
	SVec		tags;
	String		vl;
	char		*p, *var1, *var2;
	int		ncol, i, nch;

	h_hnocol( id, &ncol );

	tab = HashInt_new( ncol );
	if ( tab == 0 ) {
		return R_ALLOC_ERROR;
	}

	all_tags = svec_new( ncol );
	tags = svec_new( ncol );

	if ( HNTNEW(id) != 0 ) {
		for ( i=1 ; i <= ncol ; i++ ) {
			int	nsub, ityp, isize, ielem;
			char	tag[MAX_NAME_LEN+1], block[MAX_BLOCK_LEN+1];

			tag[0] = '\0';
			block[0] = '\0';
			HNTVAR( id, i, tag, block, nsub, ityp, isize, ielem );
			p = str_new( tag );
			svec_add( all_tags, p );
			HashInt_add( tab, p, i-1 );
		}
	} else {
		h_rwn_getInfo( id );
		for ( i=1 ; i <= ncol ; i++ ) {
			p = str_new( h_rwn_tags[i-1] );
			svec_add( all_tags, p );
			HashInt_add( tab, p, i-1 );
		}
	}

	vl = str_new( var_list );

	p = ku_getl();
	while (p != NULL) {

	  var1 = p;
	  while( (*p!='\0') && ( isalnum(*p) || (*p=='_') || (*p=='$') ) )
	    p++;
	  if ( *p == '\0' ) {
	    /* single variable */
	    svec_add( tags, str_new( var1 ) );
	  } else if ( *p == ':' ) {
	    *p = '\0';
	    p += 1;
	    var2 = p;
	    while( (*p!='\0') &&
		   ( isalnum(*p) || (*p=='_') || (*p=='$') ) )
	      p++;
	    if ( *p == '\0' ) {
	      /* var1:var2 */
	      int	err = R_NOERR;
	      
	      add_range( tab, all_tags, var1, var2,
			 tags, &err );
	      
	      if ( err != R_NOERR ) {
		svec_del( tags );
		svec_del( all_tags );
		HashInt_del( tab );
		str_del( vl );
		return R_SYNTAX_ERROR;
	      }
	    } else {
	      /* syntax error */
	      sf_report( "Illegal element in argument VARLIS"
			 " (%s:%s)\n", var1, var2 );
	      sf_report( "Should be \"expr\" or \"var1:var2\"\n");
	      svec_del( tags );
	      svec_del( all_tags );
	      HashInt_del( tab );
	      str_del( vl );
	      return R_SYNTAX_ERROR;
	    }
	  } else {
	    /* single expression */
	    svec_add( tags, str_new( var1 ) );
	  }
	p = ku_getl();
	}
	
	if ( svec_entries( tags ) == 0 ) {
		svec_del( tags );
		tags = svec_copy( all_tags );
	}

	*vlp = tags;
	svec_del( all_tags );
	HashInt_del( tab );
	str_del( vl );

	return R_NOERR;
}


/*
 *      S   Graphical scan (spider plot).
 *     ' '  Alphanumeric output of the Ntuple.
 *      S2  Graphical scan (segments plot).
 *      A   Used with 'S' it displays the average spider.
 */

static void
do_nt_scan()
{
	char		*id_string, *id_path;
	char		*selection_string;
#ifdef LONG_LONG_CHAIN
	long long int	nevent, ifirst;
#else
	int		nevent, ifirst;
#endif
	char		*opt_string;
	char		*var_lis, *tmp_lis;
	SVec		expressions;
	int		n, i, err, id;
	QuerySrc *	qs;
	QueryExe *	qe;
	bool		use_pawpp, use_spider, use_average;
	int		spider_type;

	/* get all the parameters from kuip */
	id_string = str_new( ku_getf() );
	selection_string = str_new( ku_getf() );
	C2FCBSTR( selection_string, PCCSE2.chcsel, 0 );
	nevent = ku_geti();
        if (nevent < 0) {
		/* catch old motif implementation ... */
		qp_abort( "do_nt_scan: nevent < 0 ???\n" );
	}
	ifirst = ku_geti();
	opt_string = str_tolower( ku_getc() );
	var_lis = str_new( ku_getf() );

	/* load the ntuple into memory */

	err = h_load_nt( id_string, &id_path, &id );
	str_del( id_string );
	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( selection_string );
		str_del( opt_string );
		str_del( var_lis );
		return ;
	}

	err = expand_var_list( id, var_lis, &expressions );
	str_del( var_lis );
	if ( err != R_NOERR ) {
	        h_reset_dir();
		str_del( selection_string );
		str_del( opt_string );
		return ;
	}

	qs = qp_qs_new( id_path, id, selection_string, expressions );
	str_del( id_path );
	str_del( selection_string );
	svec_del( expressions );

	qe = qp_compile( qs, FALSE, &err );
	qp_qs_free( qs );

	if ( err != R_NOERR ) {
	        h_reset_dir();
		return;
	}

	if ( ifirst < 0 ) {
		ifirst = -ifirst;
		if ( PCADDR.jmlab != 0 ) {
			use_pawpp = TRUE;
		} else {
			use_pawpp = FALSE;
		}
	} else {
		use_pawpp = FALSE;
	}

	calc_event_range( qe->id, &ifirst, &nevent );

	/* process the options */
	use_spider = FALSE;
	use_average = FALSE;

	if ( strchr( opt_string, 's' ) != 0 ) {
		use_spider = TRUE;
		spider_type = 1;
		if ( strchr( opt_string, '2' ) != 0 ) {
			spider_type = 2;
		}
	}
	if ( strchr( opt_string, 'a' ) != 0 ) {
		use_average = TRUE;
	}
	str_del( opt_string );
	
	if (use_spider && (qe->nexpr<3)) {
	  printf(" ==> Spider needs at least three variables!\n");
	  printf(" ==> Ignoring option -s ...\n");
	  use_spider = FALSE;
	}

	/* execute the command */

	if ( use_spider ) {
		QPCmd     *cmd, *cmd_tmp;
		int       zonex, zoney;

#ifdef PIAF_PROGRESS_METER
		PCNTWK.npass = 2;
		PCNTWK.ipass = 1;
		PCNTWK.misbyt = 0;
#endif
		cmd = qpcmd_new( CMD_SPIDERSCAN );
		if (!use_average) {
		  cmd_tmp=qpcmd_new( CMD_MINMAX );

		  qp_execute( qe, ifirst, nevent, cmd_tmp, &err );
		  for (i=0; i<qe->nexpr; i++) {
		    qpcmd_getminmax(cmd_tmp, i, &(cmd->u.sp_scan.min[i]), &(cmd->u.sp_scan.max[i]) );
		  }
		}
		else {
		  cmd_tmp=qpcmd_new( CMD_MINMAXAVG );

		  qp_execute( qe, ifirst, nevent, cmd_tmp, &err );
		  for (i=0; i<qe->nexpr; i++) {
		    qpcmd_getminmaxavg(cmd_tmp, i, &(cmd->u.sp_scan.min[i]), &(cmd->u.sp_scan.max[i]), &(cmd->u.sp_scan.avg[i]));
		  }
		}
		qpcmd_free( cmd_tmp );
		
		HPLGZO(zonex,zoney);
		cmd->u.sp_scan.max_line = zonex*zoney;
		cmd->u.sp_scan.rzone = (float) (zonex >= zoney ? zonex : zoney);
		cmd->u.sp_scan.spider_type = spider_type;
		cmd->u.sp_scan.use_average = use_average;
		cmd->u.sp_scan.cur_line = 0;
#ifdef PIAF_PROGRESS_METER
		PCNTWK.ipass = 2;
#endif
		qp_execute( qe, ifirst, nevent, cmd, &err );

		qpcmd_free( cmd );
	} else {
		QPCmd *		cmd;

		cmd = qpcmd_new( CMD_SCAN );
		cmd->u.scan.pawpp = use_pawpp;
		cmd->u.scan.max_line = 19; /* should be setable somewhere */
		cmd->u.scan.cur_line = 0;

#ifdef PIAF_PROGRESS_METER
		PCNTWK.npass = 1;
		PCNTWK.ipass = 1;
#endif
		qp_execute( qe, ifirst, nevent, cmd, &err );

		qpcmd_free( cmd );
	}

        h_reset_dir();

	qp_qe_free( qe );
}


static void
do_hplot_1d (
	QueryExe *	qe,
#ifdef LONG_LONG_CHAIN
	long long int	ifirst,
	long long int	nevent,
#else
	int		ifirst,
	int		nevent,
#endif
	int		idh,
	PlotOptions *	opt,
	int		*errp
){
	float		min, max;
	bool		labels = FALSE;
	SMap		label_list;

	if ( qe->expr_type[0] == D_STR ) {
		labels = TRUE;
	}

	label_list = 0;

	if ( HEXIST(idh) ) {
		if ( idh == 1000000 ) {
			HDELET( idh );
		} else {
			if ( !h_flag_1d(idh) || h_flag_profile(idh) ) {
				sf_report( "Histogram %d is not 1 dimensional,"
					" 2 expressions are required instead"
					" of 1\n", idh);
				*errp = R_SYNTAX_ERROR;
				return;
			}
			if ( labels ) {
				if ( ! HLABEQ( idh, "X" ) ) {
					sf_report( "Histogram %d has no "
						"labels\n", idh);
					*errp = R_SYNTAX_ERROR;
					return;
				}
				label_list = h_get_labels( idh, "X" );
			}
			HRESET( idh, "" );
		}
	}

#ifdef PIAF_PROGRESS_METER
	PCNTWK.npass = 1;
	PCNTWK.ipass = 0;
	PCNTWK.misbyt = 0;
#endif
	if ( ! HEXIST(idh) ) {	/* determine binning ... */

		int		nbin;
		float		min, max, low, up;
		bool		is_int;
		QPCmd *		cmd = 0;

		if (opt->S) {
			if ( labels ) {
				sf_report( "Option S not supported for "
					"character expressions\n" );
				*errp = R_SYNTAX_ERROR;
			} else {
				float *rval;
                                rval = (float *) calloc( sizeof(float), 4 ); qp_assert( rval );
				IGQWK(1,"NTWN",rval);
				if ( theLastZone.idim == 0 ) {
					is_int = ( (qe->expr_type[0] == D_INT)  || 
                        	                   (qe->expr_type[0] == D_LONG) ||
                                	           (qe->expr_type[0] == D_UINT) ||
                                        	   (qe->expr_type[0] == D_ULONG) );
					calc_axis( 100,  rval[0], rval[1], is_int,
						&low, &up, &nbin );
				} else {
					if ( theLastZone.xlow  != rval[0] ||
					     theLastZone.xhigh != rval[1]) {
						is_int = ( (qe->expr_type[0] == D_INT)  || 
                                        		   (qe->expr_type[0] == D_LONG) ||
                                    	        	   (qe->expr_type[0] == D_UINT) ||
                                                           (qe->expr_type[0] == D_ULONG) );
						calc_axis( 100,  rval[0], rval[1], is_int,
							  &low, &up, &nbin );
					} else {
						nbin = theLastZone.nxbin;
						low  = theLastZone.xlow;
						up   = theLastZone.xhigh;
					}
				}
				free( (void *) rval );
			}
		} else {

			cmd = qpcmd_new( CMD_MINMAX );

#ifdef PIAF_PROGRESS_METER
			PCNTWK.npass = 2;
			PCNTWK.ipass = 1;
#endif
			qp_execute( qe, ifirst, nevent, cmd, errp );

			if ((HCPIAF.slavpf) && (*errp != R_NOERR)) {
			  PiafEmpty1D (qe->expr_str[0],labels);
			}
			  
			if ((!cmd->u.minmax.has_events)&&(*errp == R_NOERR)){
				*errp = R_NOEVT;
				if ( HCPIAF.slavpf ) {
				  PiafEmpty1D (qe->expr_str[0],labels);
				} else {
				  if (!opt->N) {
					HPLFRA(0.,10.,0.,10.,"A");
					IGTEXT(5.,5.,"Empty",.3,text_angle,"C");
					HPLTIT(qe->expr_str[0]);
				  }
				  else {
				    sf_report("Empty Histogram not created\n");
				  }
				}
			}

			if ( *errp == R_NOERR && !labels ) {

				qpcmd_getminmax( cmd, 0, &min, &max );
				if ( HCPIAF.slavpf ) {
				  PFMINMAX (min,max);
				}

				is_int = ( (qe->expr_type[0] == D_INT) || 
					(qe->expr_type[0] == D_LONG) ||
					(qe->expr_type[0] == D_UINT) ||
					(qe->expr_type[0] == D_ULONG) );

				calc_axis( 100, min, max, is_int,
					&low, &up, &nbin );

				theLastZone.idim = 1;
				theLastZone.nxbin = nbin;
				theLastZone.xlow = low;
				theLastZone.xhigh = up;
			}

			if ( *errp == R_NOERR && labels ) {

				label_list = qpcmd_labels( cmd, 0 );
				if ( HCPIAF.slavpf ) {
					String s;
					int len = 32768, pos = 0;
					s = str_alloc( len + 1 );
					PackLabels(label_list,&pos,s,&len);
					PFLABELS( s, len);
					UnpackLabels( label_list, s, len );
					str_del( s );
				}
				theLastZone.idim = 0;
				theLastZone.nxbin = 0;
				theLastZone.xlow = 0;
				theLastZone.xhigh = 0;
			}
			qpcmd_free( cmd );
		}

		if ( *errp == R_NOERR ) {
			if ( labels ) {
				h_hbook1_labels( idh, qe->expr_str[0],
					label_list );
			} else {
				HBOOK1( idh, qe->expr_str[0],
					nbin, low, up, 0. );
			}
		}

	}

	if ( *errp == R_NOERR ) {
		QPCmd *		cmd;
		QPCmdHFill1	*h;
		
		cmd = qpcmd_new( CMD_HFILL1 );
		h = &cmd->u.hfill1;

		h->idh = idh;
		h->cvt_x = datatype_to_cvtcallback( qe->expr_type[0],
				(void *) label_list );

#ifdef PIAF_PROGRESS_METER
		PCNTWK.ipass++;
#endif
		qp_execute( qe, ifirst, nevent, cmd, errp );

		qpcmd_free( cmd );
	}

	if (( *errp == R_NOERR ) && ( !opt->N )) {
		char *		opt_string;

		opt_string = qp_plot_opt_gen( opt, TRUE );
		PAPLOT( idh, opt_string, "", 0, 0, 0, 0, 0, 0 );
		PAUTIT( " " );

		str_del( opt_string );
	}

	if ( label_list != 0 ) {
		smap_del( label_list );
	}
}


static void
do_hplot_2d (
	QueryExe *	qe,
#ifdef LONG_LONG_CHAIN
	long long int	ifirst,
	long long int	nevent,
#else
	int		ifirst,
	int		nevent,
#endif
	int		idh,
	PlotOptions *	opt,
	int *		errp
) {
	bool		labelx = FALSE, labely = FALSE;
	SMap		labelx_list, labely_list;

	labelx_list = 0;
	labely_list = 0;

	if ( qe->expr_type[0] == D_STR ) {
		labely = TRUE;
	}

	if ( qe->expr_type[1] == D_STR ) {
		labelx = TRUE;
	}

	if ( HEXIST(idh) ) {
		if ( idh == 1000000 ) {
			HDELET( idh );
		} else {
			if ( !h_flag_2d(idh) && !h_flag_profile(idh) ) {
				sf_report( "Histogram %d is not 2-dimensional"
					"or profile\n", idh);
				*errp = R_SYNTAX_ERROR;
				return;
			}

			if ( labelx ) {
				if ( ! HLABEQ( idh, "X" ) ) {
					sf_report( "Histogram %d has no "
						"labels on X axis\n", idh);
					*errp = R_SYNTAX_ERROR;
					return;
				}
				labelx_list = h_get_labels( idh, "X" );
			}

			if ( labely ) {
				if ( ! HLABEQ( idh, "Y" ) ) {
					sf_report( "Histogram %d has no "
						"labels on Y axis\n", idh);
					*errp = R_SYNTAX_ERROR;
					return;
				}
				labely_list = h_get_labels( idh, "Y" );
			}

			HRESET( idh, "" );
		}
	}

#ifdef PIAF_PROGRESS_METER
	PCNTWK.npass = 1;
	PCNTWK.ipass = 0;
	PCNTWK.misbyt = 0;
#endif

	if ( ! HEXIST(idh) ) {	/* determine binning ... */

		int		nbin[2], i;
		float		min, max, low[2], up[2];
		bool		is_int;
		char *		title;
		QPCmd *		cmd;

		title = str_merge(
			str_new(qe->expr_str[0]),
			str_new(" VS. "),
			str_new(qe->expr_str[1]),
			(char *) 0 );

		if ( opt->S ) {
			if ( labelx || labely ) {
				sf_report( "Option S not supported for "
					"character expressions\n" );
				*errp = R_SYNTAX_ERROR;
			} else if ( theLastZone.idim == 0 ) {
				sf_report( "No info about previous plot\n" );
				*errp = R_SYNTAX_ERROR;
			} else if ( theLastZone.idim != 2 ) {
				sf_report( "Cannot overlay a 2d histogram on a"
					" 1d zone\n" );
				*errp = R_SYNTAX_ERROR;
			} else {
				nbin[1] = theLastZone.nxbin;
				low[1] = theLastZone.xlow;
				up[1] = theLastZone.xhigh;
				nbin[0] = theLastZone.nybin;
				low[0] = theLastZone.ylow;
				up[0] = theLastZone.yhigh;
			}
		} else {

			cmd = qpcmd_new( CMD_MINMAX );

#ifdef PIAF_PROGRESS_METER
			PCNTWK.npass = 2;
			PCNTWK.ipass = 1;
#endif
			qp_execute( qe, ifirst, nevent, cmd, errp );

			if ((HCPIAF.slavpf) && (*errp != R_NOERR)) {
			  PiafEmpty2D(title,labelx,labely);
			}

			if ((!cmd->u.minmax.has_events)&&(*errp == R_NOERR)){
				*errp = R_NOEVT;
				if ( HCPIAF.slavpf ) {
					PiafEmpty2D(title,labelx,labely);
				} else {
				  if (!opt->N) {
					HPLFRA(0.,10.,0.,10.,"A");
					IGTEXT(5.,5.,"Empty",.3,text_angle,"C");
					HPLTIT( title );
				  }
				  else {
				    sf_report("Empty Histogram not created\n");
				  }
				}
			}

			if ( *errp == R_NOERR ) {

				for( i=0 ; i < 2 ; i++ ) {

					if ( i==0 && labely ) {
						labely_list = qpcmd_labels(
								cmd, 0 );
						if ( HCPIAF.slavpf ) {
							String s;
							int len = 32768;
							int pos = 0;
							s = str_alloc(len+1);
							PackLabels(
								labely_list, 
								&pos, s,&len);
							PFLABELS( s, len);
							UnpackLabels( 
								labely_list,
								s, len );
							str_del( s );
						}

						continue;
					}

					if ( i==1 && labelx ) {
						labelx_list = qpcmd_labels(
								cmd, 1 );
						if ( HCPIAF.slavpf ) {
							String s;
							int len = 32768;
							int pos = 0;
							s = str_alloc(len+1);
							PackLabels(
								labelx_list, 
								&pos, s,&len);
							PFLABELS( s, len);
							UnpackLabels( 
								labelx_list,
								s, len );
							str_del( s );
						}

						continue;
					}

					qpcmd_getminmax( cmd, i, &min, &max );

					if ( HCPIAF.slavpf ) {
					  PFMINMAX (min,max);
					}

					is_int =(
						(qe->expr_type[i] == D_INT) || 
						(qe->expr_type[i] == D_LONG) ||
						(qe->expr_type[i] == D_UINT) ||
						(qe->expr_type[i] == D_ULONG) );

					calc_axis( 40, min, max, is_int,
						&low[i], &up[i], &nbin[i] );
				}

				if ( labelx || labely ) {
					theLastZone.idim = 0;
				} else {
					theLastZone.idim = 2;
					theLastZone.nxbin = nbin[1];
					theLastZone.xlow = low[1];
					theLastZone.xhigh = up[1];
					theLastZone.nybin = nbin[0];
					theLastZone.ylow = low[0];
					theLastZone.yhigh = up[0];
				}
			}

			qpcmd_free( cmd );
		}

		if ( *errp == R_NOERR ) {

			h_hbook2_labels( idh, title,
				labelx_list, labely_list,
				nbin, low, up );
			
		}

		str_del( title );
	}


	if ( *errp == R_NOERR ) {
		QPCmd *		cmd;
		QPCmdHFill2	*h;

		cmd = qpcmd_new( CMD_HFILL2 );
		h = &cmd->u.hfill2;


		h->idh = idh;
		h->cvt_x = datatype_to_cvtcallback( qe->expr_type[1],
							labelx_list );
		h->cvt_y = datatype_to_cvtcallback( qe->expr_type[0],
							labely_list );

#ifdef PIAF_PROGRESS_METER
		PCNTWK.ipass++;
#endif
		qp_execute( qe, ifirst, nevent, cmd, errp );

		qpcmd_free( cmd );
	}


	if (( *errp == R_NOERR ) && ( !opt->N )) {
		char *		opt_string;

		opt_string = qp_plot_opt_gen( opt, TRUE );
		PAPLOT( idh, opt_string, "", 0, 0, 0, 0, 0, 0 );
		PAUTIT( " " );

		str_del( opt_string );
	}

	if ( labelx_list != 0 ) {
		smap_del( labelx_list );
	}
	if ( labely_list != 0 ) {
		smap_del( labely_list );
	}
}


static void
do_prof_1d (
	QueryExe *	qe,
#ifdef LONG_LONG_CHAIN
	long long int	ifirst,
	long long int	nevent,
#else
	int		ifirst,
	int		nevent,
#endif
	int		idh,
	PlotOptions *	opt,
	int		*errp
) {

	if ( (qe->expr_type[0] == D_STR) || (qe->expr_type[1] == D_STR) ) {
		sf_report( "do_prof_1d: D_STR not implemented\n" );
		*errp = R_NOT_IMPLEMENTED;
		return;
	}

	if ( HEXIST(idh) ) {
		if ( idh == 1000000 ) {
			HDELET( idh );
		} else if ( h_flag_profile(idh) ) {
			HRESET( idh, "" );
		} else {
			sf_report( "Histogram %d is not a profile histogram\n", idh );
			*errp = R_TYPE_ERROR;
			return;
		}
	}

#ifdef PIAF_PROGRESS_METER
	PCNTWK.npass = 1;
	PCNTWK.ipass = 0;
	PCNTWK.misbyt = 0;
#endif
	if ( ! HEXIST(idh) ) {	/* determine binning ... */

		int		nbin[2], i;
		float		min, max, low[2], up[2];
		bool		is_int;
		char *		title;
		QPCmd *		cmd;

		title = str_merge(
			str_new(qe->expr_str[0]),
			str_new(" VS. "),
			str_new(qe->expr_str[1]),
			(char *) 0 );

		if ( opt->S ) {
			if ( theLastZone.idim != 2 ) {
				sf_report( "Cannot overlay a 2d histogram on a"
					" 1d zone\n" );
				*errp = R_SYNTAX_ERROR;
			} else {
				nbin[1] = theLastZone.nxbin;
				low[1] = theLastZone.xlow;
				up[1] = theLastZone.xhigh;
				nbin[0] = theLastZone.nybin;
				low[0] = theLastZone.ylow;
				up[0] = theLastZone.yhigh;
			}
		} else {
			cmd = qpcmd_new( CMD_MINMAX );

#ifdef PIAF_PROGRESS_METER
			PCNTWK.npass = 2;
			PCNTWK.ipass = 1;
#endif
			qp_execute( qe, ifirst, nevent, cmd, errp );

			if ((HCPIAF.slavpf) && (*errp != R_NOERR)) {
			  PiafEmptyNoLabels(title,2);
			}

			if ((!cmd->u.minmax.has_events)&&(*errp == R_NOERR)){
				*errp = R_NOEVT;
				if ( HCPIAF.slavpf ) {
					PiafEmptyNoLabels(title,2);
				} else {
				  if (!opt->N) {
					HPLFRA(0.,10.,0.,10.,"A");
					IGTEXT(5.,5.,"Empty",.3,text_angle,"C");
					HPLTIT( title );
				  }
				  else {
				    sf_report("Empty Histogram not created\n");
				  }
				}
			}

			if ( *errp == R_NOERR ) {

				for( i=0 ; i < 2 ; i++ ) {

					qpcmd_getminmax( cmd, i, &min, &max );

					if ( HCPIAF.slavpf ) {
					  PFMINMAX (min,max);
					}

					is_int =(
						(qe->expr_type[i] == D_INT) || 
						(qe->expr_type[i] == D_LONG) ||
						(qe->expr_type[i] == D_UINT) ||
						(qe->expr_type[i] == D_ULONG) );

					calc_axis( 100, min, max, is_int,
						&low[i], &up[i], &nbin[i] );
				}

				theLastZone.idim = 2;
				theLastZone.nxbin = nbin[1];
				theLastZone.xlow = low[1];
				theLastZone.xhigh = up[1];
				theLastZone.nybin = nbin[0];
				theLastZone.ylow = low[0];
				theLastZone.yhigh = up[0];

			}

			qpcmd_free( cmd );
		}

		if ( *errp == R_NOERR ) {
			char	*optstr;

			if ( opt->profi ) {
				optstr = "I";
			} else if ( opt->profs ) {
				optstr = "S";
			} else {
				optstr = " ";
			}

			HBPROF( idh, title, nbin[1],
				low[1], up[1],
				low[0], up[0], optstr );
			
		}

		str_del( title );
	}

	if ( *errp == R_NOERR ) {
		QPCmd *		cmd;
		QPCmdHFill2	*h;

		cmd = qpcmd_new( CMD_HFILL2 );
		h = &cmd->u.hfill2;

		h->idh = idh;
		h->cvt_x = datatype_to_cvtcallback( qe->expr_type[1],
							(void *) 0 );
		h->cvt_y = datatype_to_cvtcallback( qe->expr_type[0],
							(void *) 0 );

#ifdef PIAF_PROGRESS_METER
		PCNTWK.ipass++;
#endif
		qp_execute( qe, ifirst, nevent, cmd, errp );

		qpcmd_free( cmd );
	}

	if (( *errp == R_NOERR ) && ( !opt->N )) {
		char *		opt_string;

		opt_string = qp_plot_opt_gen( opt, TRUE );
		PAPLOT( idh, opt_string, "", 0, 0, 0, 0, 0, 0 );
		PAUTIT( " " );

		str_del( opt_string );
	}
}


static void
do_scat_2d (
	QueryExe *	qe,
#ifdef LONG_LONG_CHAIN
	long long int	ifirst,
	long long int	nevent,
#else
	int		ifirst,
	int		nevent,
#endif
	int		idh,
	PlotOptions *	opt,
	int		*errp
) {
	bool		labelx = FALSE, labely = FALSE;
	SMap		labelx_list, labely_list;

	labelx_list = 0;
	labely_list = 0;

	if ( qe->expr_type[0] == D_STR ) {
		labely = TRUE;
	}

	if ( qe->expr_type[1] == D_STR ) {
		labelx = TRUE;
	}

	if ( HEXIST(idh) ) {
		if ( idh == 1000000 ) {
			HDELET( idh );
		} else {
			if ( !h_flag_2d(idh) && !h_flag_profile(idh) ) {
				sf_report( "Histogram %d is not 2-dimensional"
					"or profile\n", idh);
				*errp = R_SYNTAX_ERROR;
				return;
			}

			if ( labelx ) {
				if ( ! HLABEQ( idh, "X" ) ) {
					sf_report( "Histogram %d has no "
						"labels on X axis\n", idh);
					*errp = R_SYNTAX_ERROR;
					return;
				}
				labelx_list = h_get_labels( idh, "X" );
			}

			if ( labely ) {
				if ( ! HLABEQ( idh, "Y" ) ) {
					sf_report( "Histogram %d has no "
						"labels on Y axis\n", idh);
					*errp = R_SYNTAX_ERROR;
					return;
				}
				labely_list = h_get_labels( idh, "Y" );
			}

			HRESET( idh, "" );
		}
	}

#ifdef PIAF_PROGRESS_METER
	PCNTWK.npass = 1;
	PCNTWK.ipass = 0;
	PCNTWK.misbyt = 0;
#endif
	if ( ! HEXIST(idh) ) {	/* determine binning ... */

		int		nbin[2], i;
		float		min, max, low[2], up[2];
		bool		is_int;
		char *		title;
		QPCmd *		cmd;

		title = str_merge(
			str_new(qe->expr_str[0]),
			str_new(" VS. "),
			str_new(qe->expr_str[1]),
			(char *) 0 );

		if ( opt->S ) {
			if ( labelx || labely ) {
				sf_report( "Option S not supported for "
					"character expressions\n" );
				*errp = R_SYNTAX_ERROR;
			} else {
				float *rval;
				rval = (float *) calloc( sizeof(float), 4 ); qp_assert( rval );
				IGQWK(1,"NTWN",rval);
				if ( theLastZone.idim == 0 ) {
					nbin[1] = 40;
					low[1]  = rval[0];
					up[1]   = rval[1];
					nbin[0] = 40;
					low[0]  = rval[2];
					up[0]   = rval[3];
				} else {
					if ( theLastZone.xlow  != rval[0] ||
					     theLastZone.xhigh != rval[1] ||
					     theLastZone.ylow  != rval[2] ||
					     theLastZone.yhigh != rval[3] ) {
						nbin[1] = 40;
						low[1]  = rval[0];
						up[1]   = rval[1];
						nbin[0] = 40;
						low[0]  = rval[2];
						up[0]   = rval[3];
					} else {
						nbin[1] = theLastZone.nxbin;
						low[1]  = theLastZone.xlow;
						up[1]   = theLastZone.xhigh;
						nbin[0] = theLastZone.nybin;
						low[0]  = theLastZone.ylow;
						up[0]   = theLastZone.yhigh;
					}
				}
				free( (void *) rval );
			}
		} else {
			cmd = qpcmd_new( CMD_MINMAX );

#ifdef PIAF_PROGRESS_METER
			PCNTWK.npass = 2;
			PCNTWK.ipass = 1;
#endif
			qp_execute( qe, ifirst, nevent, cmd, errp );

			if ((HCPIAF.slavpf) && (*errp != R_NOERR)) {
			  PiafEmpty2D(title,labelx,labely);
			}

			if ((!cmd->u.minmax.has_events)&&(*errp == R_NOERR)){
				*errp = R_NOEVT;
				if ( HCPIAF.slavpf ) {
					PiafEmpty2D(title,labelx,labely);
				} else {
					HPLFRA(0.,10.,0.,10.,"A");
					IGTEXT(5.,5.,"Empty",.3,text_angle,"C");
					HPLTIT( title );
				}
			}

			if ( *errp == R_NOERR ) {

				for( i=0 ; i < 2 ; i++ ) {

					if ( i==0 && labely ) {
						labely_list = qpcmd_labels(
								cmd, 0 );
						if ( HCPIAF.slavpf ) {
							String s;
							int len = 32768;
							int pos = 0;
							s = str_alloc(len+1);
							PackLabels(
								labely_list, 
								&pos, s,&len);
							PFLABELS( s, len);
							UnpackLabels( 
								labely_list,
								s, len );
							str_del( s );
						}

						continue;
					}

					if ( i==1 && labelx ) {
						labelx_list = qpcmd_labels(
								cmd, 1 );
						if ( HCPIAF.slavpf ) {
							String s;
							int len = 32768;
							int pos = 0;
							s = str_alloc(len+1);
							PackLabels(
								labelx_list, 
								&pos, s,&len);
							PFLABELS( s, len);
							UnpackLabels( 
								labelx_list,
								s, len );
							str_del( s );
						}

						continue;
					}

					qpcmd_getminmax( cmd, i, &min, &max );
			  
					if ( HCPIAF.slavpf ) {
					  PFMINMAX (min,max);
					}

					is_int =(
						(qe->expr_type[i] == D_INT) || 
						(qe->expr_type[i] == D_LONG) ||
						(qe->expr_type[i] == D_UINT) ||
						(qe->expr_type[i] == D_ULONG) );

					calc_axis( 40, min, max, is_int,
						&low[i], &up[i], &nbin[i] );
				}

				if ( labelx || labely ) {
					theLastZone.idim = 0;
				} else {
					theLastZone.idim = 2;
					theLastZone.nxbin = nbin[1];
					theLastZone.xlow = low[1];
					theLastZone.xhigh = up[1];
					theLastZone.nybin = nbin[0];
					theLastZone.ylow = low[0];
					theLastZone.yhigh = up[0];
				}
			}

			qpcmd_free( cmd );
		}

		if ( *errp == R_NOERR ) {

			h_hbook2_labels( idh, title,
				labelx_list, labely_list,
				nbin, low, up );
		}
			
		str_del( title );
	}

	if ( *errp == R_NOERR ) {
		QPCmd *		cmd;
		QPCmdScat2	*h;
                bool		stats_required;
		char *		opt_string;

		opt_string = qp_plot_opt_gen( opt, TRUE );
		cmd = qpcmd_new( CMD_SCAT2 );
		h = &cmd->u.scat2;

		h->idh = idh;
		h->line = opt->line;
		h->logx = opt->logx;
		h->logy = opt->logy;

		h->cvt_x = datatype_to_cvtcallback( qe->expr_type[1],
							labelx_list );
		h->cvt_y = datatype_to_cvtcallback( qe->expr_type[0],
							labely_list );

		HPLOPT("STA ",-1);
		stats_required = QUEST.iquest[10] != 0;

		if ( stats_required ) {
			HPLOPT( "NSTA", 1 );
		}

                PAPLOT( idh, opt_string, " ", 0, 0, 0, 0, 0, 0 );
		if ( HCPIAF.slavpf ) {
			char buf[80];
			int istat;
			strcpy (buf,"PSYNC");
			CZPUTA (buf, istat);
			CZGETA (buf, istat);
		}

		PAUTIT( " " );
                IGPID(1,"ntuple",qe->id," ");

#ifdef PIAF_PROGRESS_METER
		PCNTWK.ipass++;
#endif
		qp_execute( qe, ifirst, nevent, cmd, errp );

	        if ( stats_required ) {
			HFIND(idh,"do_scat_2d");
			HDCOFL();
			HPLOPT( "STA ", 1 );
			HPLSTA( idh, "HIST", 1 );
		}

		str_del( opt_string );
		qpcmd_free( cmd );
	}

	if ( labelx_list != 0 ) {
		smap_del( labelx_list );
	}
	if ( labely_list != 0 ) {
		smap_del( labely_list );
	}
}


static void
do_scat_3d (
	QueryExe *	qe,
#ifdef LONG_LONG_CHAIN
	long long int	ifirst,
	long long int	nevent,
#else
	int		ifirst,
	int		nevent,
#endif
	PlotOptions *	opt,
	int		*errp
) {
	int		nbin[3], i;
	float		min, max, low[3], up[3], Theta, Phi;
	bool		is_int;
	char *		title;
	QPCmd *		cmd;

/*
 *  Scatter plot on character variables not yet implemented
 */

	if ( (qe->expr_type[0] == D_STR) ||
		(qe->expr_type[1] == D_STR) || (qe->expr_type[2] == D_STR) ) {
		sf_report( "Character expressions not available for 3 "
			"dimensional plots\n" );
		*errp = R_NOT_IMPLEMENTED;
		return;
	}

	title = str_merge(
		str_new(qe->expr_str[0]),
		str_new(" VS. "),
		str_new(qe->expr_str[1]),
		str_new(" VS. "),
		str_new(qe->expr_str[2]),
		(char *) 0 );

#ifdef PIAF_PROGRESS_METER
	PCNTWK.npass = 1;
	PCNTWK.ipass = 0;
	PCNTWK.misbyt = 0;
#endif
	if(!opt->S) {
	  cmd = qpcmd_new( CMD_MINMAX );
	  
#ifdef PIAF_PROGRESS_METER
	  PCNTWK.npass = 2;
	  PCNTWK.ipass = 1;
#endif
	  qp_execute( qe, ifirst, nevent, cmd, errp );
	  
	  if ((HCPIAF.slavpf) && (*errp != R_NOERR)) {
		PiafEmptyNoLabels(title,3);
	  }

	  if ((!cmd->u.minmax.has_events)&&(*errp == R_NOERR)){
	    *errp = R_NOEVT;
	    if ( HCPIAF.slavpf ) {
	      PiafEmptyNoLabels(title,3);
	    } else {
	      HPLFRA(0.,10.,0.,10.,"A");
	      IGTEXT(5.,5.,"Empty",.3,text_angle,"C");
	      HPLTIT( title );
	    }
	  }
	  
	  if ( *errp == R_NOERR ) {
	    
	    for( i=0 ; i < 3 ; i++ ) {
	      
	      qpcmd_getminmax( cmd, i, &min, &max );
	      
	      if ( HCPIAF.slavpf ) {
		PFMINMAX (min,max);
	      }
	      
	      is_int = ( (qe->expr_type[i] == D_INT) || 
			 (qe->expr_type[i] == D_LONG) ||
			 (qe->expr_type[i] == D_UINT) ||
			 (qe->expr_type[i] == D_ULONG) );
	      
	      calc_axis( 40, min, max, is_int,
			 &low[i], &up[i], &nbin[i] );
	    }
	    
	    Theta = 0.;
	    Phi   = 0.;
	    GETTP( Theta, Phi );
	    HPLFR3( low[2], up[2], low[1], up[1], low[0], up[0],
		    Theta, Phi, "BW" );
	    
	  }
	  qpcmd_free( cmd );
	}

	PAUTIT( title );
	str_del( title );

	if ( *errp == R_NOERR ) {
		QPCmd *		cmd;

		cmd = qpcmd_new( CMD_SCAT3 );
		cmd->u.scat3.line = opt->line;
		cmd->u.scat3.logx = opt->logx;
		cmd->u.scat3.logy = opt->logy;
		cmd->u.scat3.logz = opt->logz;

                IGPID(1,"ntuple",qe->id," ");

#ifdef PIAF_PROGRESS_METER
		PCNTWK.ipass++;
#endif
		qp_execute( qe, ifirst, nevent, cmd, errp );

		if (!opt->S) {
		  HPLFR3( low[2], up[2], low[1], up[1], low[0], up[0],
			Theta, Phi, "F" );
		}

		qpcmd_free( cmd );
	}

}

static void
do_gouraud (
	QueryExe *	qe,
#ifdef LONG_LONG_CHAIN
	long long int	ifirst,
	long long int	nevent,
#else
	int		ifirst,
	int		nevent,
#endif
	PlotOptions *	opt,
	int		*errp
) {
	int		nbin[3], i;
	float		min, max, low[3], up[3], Theta, Phi;
	bool		is_int;
	char *		title;
	QPCmd *		cmd_tmp;
	QPCmd *		cmd;

/*
 *  Gouraud  plot on character variables not yet implemented
 */

	if ( (qe->expr_type[0] == D_STR) ||
		(qe->expr_type[1] == D_STR) || (qe->expr_type[2] == D_STR) ) {
		sf_report( "Character expressions not available for 3 "
			"dimensional plots\n" );
		*errp = R_NOT_IMPLEMENTED;
		return;
	}

	title = str_merge(
		str_new(qe->expr_str[0]),
		str_new(" VS. "),
		str_new(qe->expr_str[1]),
		str_new(" VS. "),
		str_new(qe->expr_str[2]),
		(char *) 0 );

	cmd_tmp = qpcmd_new( CMD_MINMAX );

#ifdef PIAF_PROGRESS_METER
	PCNTWK.npass = 2;
	PCNTWK.ipass = 1;
	PCNTWK.misbyt = 0;
#endif
	qp_execute( qe, ifirst, nevent, cmd_tmp, errp );

	if ((HCPIAF.slavpf) && (*errp != R_NOERR)) {
		PiafEmptyNoLabels(title,3);
	}

	if ((!cmd_tmp->u.minmax.has_events)&&(*errp == R_NOERR)){
		*errp = R_NOEVT;
		if ( HCPIAF.slavpf ) {
			PiafEmptyNoLabels(title,3);
		} else {
			HPLFRA(0.,10.,0.,10.,"A");
			IGTEXT(5.,5.,"Empty",.3,text_angle,"C");
			HPLTIT( title );
		}
	}

	if ( *errp == R_NOERR ) {

	  cmd = qpcmd_new( CMD_GOURAUD );
	  for( i=0 ; i < 3 ; i++ ) {

			qpcmd_getminmax( cmd_tmp, i, &(cmd->u.gouraud.min[2-i]), &(cmd->u.gouraud.max[2-i]) );

			if ( HCPIAF.slavpf ) {
			  PFMINMAX (cmd->u.gouraud.min[2-i], cmd->u.gouraud.max[2-i]);
			}

			is_int = ( (qe->expr_type[i] == D_INT) || 
				(qe->expr_type[i] == D_LONG) ||
				(qe->expr_type[i] == D_UINT) ||
				(qe->expr_type[i] == D_ULONG) );

				calc_axis( 40, cmd->u.gouraud.min[2-i], cmd->u.gouraud.max[2-i], is_int,
				&low[i], &up[i], &nbin[i] );
				cmd->u.gouraud.min[2-i] = low[i];
				cmd->u.gouraud.max[2-i] = up[i];
		}

		Theta = 0.;
		Phi   = 0.;
		GETTP( Theta, Phi );
		HPLFR3( low[2], up[2], low[1], up[1], low[0], up[0],
			Theta, Phi, "BW" );
		PAUTIT( title );
	}

	str_del( title );

	if ( *errp == R_NOERR ) {

		cmd->u.gouraud.logx = opt->logx;
		cmd->u.gouraud.logy = opt->logy;
		cmd->u.gouraud.logz = opt->logz;

                IGPID(1,"ntuple",qe->id," ");

#ifdef PIAF_PROGRESS_METER
		PCNTWK.ipass++;
#endif
		qp_execute( qe, ifirst, nevent, cmd, errp );
		HPLFR3( low[2], up[2], low[1], up[1], low[0], up[0],
			Theta, Phi, "FG" );

		qpcmd_free( cmd );
	}

	qpcmd_free( cmd_tmp );
}

static void
do_scat_4d (
	QueryExe *	qe,
#ifdef LONG_LONG_CHAIN
	long long int	ifirst,
	long long int	nevent,
#else
	int		ifirst,
	int		nevent,
#endif
	PlotOptions *	opt,
	int		*errp
) {
	int		nbin[4], i;
	float		min, max, low[4], up[4], Theta, Phi;
	bool		is_int;
	char *		title;
	QPCmd *		cmd;

/*
 *  Scatter plot on character variables not yet implemented
 */

	if ( (qe->expr_type[0] == D_STR) || (qe->expr_type[1] == D_STR) ||
		(qe->expr_type[2] == D_STR) || (qe->expr_type[3] == D_STR) ) {
		sf_report( "Character expressions not available for 4 "
			"dimensional plots\n" );
		*errp = R_NOT_IMPLEMENTED;
		return;
	}

	title = str_merge(
		str_new(qe->expr_str[0]),
		str_new(" VS. "),
		str_new(qe->expr_str[1]),
		str_new(" VS. "),
		str_new(qe->expr_str[2]),
		str_new(" VS. "),
		str_new(qe->expr_str[3]),
		(char *) 0 );

	cmd = qpcmd_new( CMD_MINMAX );
	
#ifdef PIAF_PROGRESS_METER
	PCNTWK.npass = 2;
	PCNTWK.ipass = 1;
	PCNTWK.misbyt = 0;
#endif
	qp_execute( qe, ifirst, nevent, cmd, errp );
	
	if ((HCPIAF.slavpf) && (*errp != R_NOERR)) {
		PiafEmptyNoLabels(title,4);
	}

	if ((!cmd->u.minmax.has_events)&&(*errp == R_NOERR)){
	  *errp = R_NOEVT;
	  if ( HCPIAF.slavpf ) {
	    PiafEmptyNoLabels(title,4);
	  } else {
	    if (!opt->S) HPLFRA(0.,10.,0.,10.,"A");
	    IGTEXT(5.,5.,"Empty",.3,text_angle,"C");
	    HPLTIT( title );
	  }
	}
	  
	  if ( *errp == R_NOERR ) {
	    QPCmd *		cmd2;
	    
	    cmd2 = qpcmd_new( CMD_SCAT4 );
	    
	    if (!opt->S) {
	      for( i=0 ; i < 4 ; i++ ) {
		
		qpcmd_getminmax( cmd, i, &min, &max );
		
		if ( HCPIAF.slavpf ) {
		  PFMINMAX (min,max);
		}
		
		cmd2->u.scat4.col_min = min;
		cmd2->u.scat4.col_max = max;
		
		is_int = ( (qe->expr_type[i] == D_INT) || 
			   (qe->expr_type[i] == D_LONG) ||
			   (qe->expr_type[i] == D_UINT) ||
			   (qe->expr_type[i] == D_ULONG) );
		
		calc_axis( 40, min, max, is_int,
			   &low[i], &up[i], &nbin[i] );
	      }
	      Theta = 0.;
	      Phi   = 0.;
	      GETTP( Theta, Phi );
	      HPLFR3( low[2], up[2], low[1], up[1], low[0], up[0],
		      Theta, Phi, "BW" );
	    }
	    else {
		qpcmd_getminmax( cmd, 3, &min, &max );
		is_int = ( (qe->expr_type[3] == D_INT) || 
			   (qe->expr_type[3] == D_LONG) ||
			   (qe->expr_type[3] == D_UINT) ||
			   (qe->expr_type[3] == D_ULONG) );
		
		calc_axis( 40, min, max, is_int,
			   &low[3], &up[3], &nbin[3] );
	    }
	      
	    cmd2->u.scat4.col_min = low[3];
	    cmd2->u.scat4.col_max = up[3];
	    cmd2->u.scat4.line = opt->line;
	    cmd2->u.scat4.logx = opt->logx;
	    cmd2->u.scat4.logy = opt->logy;
	    cmd2->u.scat4.logz = opt->logz;
	    
	    PAUTIT( title );
	    
	    IGPID(1,"ntuple",qe->id," ");
	    
#ifdef PIAF_PROGRESS_METER
	    PCNTWK.ipass++;
#endif
	    qp_execute( qe, ifirst, nevent, cmd2, errp );
	    
	    if(!opt->S) {
	      HPLFR3( low[2], up[2], low[1], up[1], low[0], up[0],
		    Theta, Phi, "F" );
	    }
	    
	    qpcmd_free( cmd2 );
	  }
	  
	  str_del( title );
	  qpcmd_free( cmd );
}


static void
do_nt_loop ( void )
{
	char		*selection_string;
	char		*id_string, *id_path;
#ifdef LONG_LONG_CHAIN
	long long int	nevent, ifirst;
#else
	int		nevent, ifirst;
#endif
	int		err;
	int		id;
	SVec		sv;
	QuerySrc *	qs;
	QueryExe *	qe;

	/* get all the parameters from kuip */
	id_string = str_new( ku_getf() );
	selection_string = str_new( ku_getf() );
	C2FCBSTR( selection_string, PCCSE2.chcsel, 0 );
	nevent = ku_geti();
	ifirst = ku_geti();

	/* load the ntuple into memory */

	err = h_load_nt( id_string, &id_path, &id );

	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( id_string );
		str_del( selection_string );
		return;
	}

	sv = svec_new( 0 );
	qs = qp_qs_new( id_path, id, selection_string, sv );
	str_del( id_path );
	svec_del( sv );

	qe = qp_compile( qs, FALSE, &err );
	qp_qs_free( qs );

	if ( err != R_NOERR ) {
	        h_reset_dir();
		str_del( id_string );
		str_del( selection_string );
		return;
	}

	calc_event_range( qe->id, &ifirst, &nevent );

	if ( HCPIAF.ntpiaf == CF_TRUE ) {
		int	istat;

		istat = 0;
		PFPING( 0, 1, istat );

		if ( istat == 0 ) {
			PFPUSH( istat );
		}
		if ( istat == 0 ) {
			char	*buf;
			int	len;

			/* rebuild command string */

			len = 33 + 1;
			len += strlen( id_string );
			len += strlen( selection_string );

			buf = (char *) calloc( len, 1 ); 
			qp_assert( buf != 0 );

#ifdef LONG_LONG_CHAIN
			sprintf( buf, "nt/loop %s %s%12lld%12lld",
#else
			sprintf( buf, "nt/loop %s %s%12d%12d",
#endif
				id_string, selection_string,
				nevent, ifirst );

			if ( qp_flags_get( "verbose" ) != 0 ) {
				sf_report( "PFKUIP(%s)\n",buf );
			}
			PFKUIP( buf, istat );

			free( (void *) buf );
		}

	} else {
		QPCmd *		cmd;

		cmd = qpcmd_new( CMD_NULL );

#ifdef PIAF_PROGRESS_METER
		PCNTWK.npass = 1;
		PCNTWK.ipass = 1;
		PCNTWK.misbyt = 0;
#endif
		qp_execute( qe, ifirst, nevent, cmd, &err );

		qpcmd_free( cmd );
	}

        h_reset_dir();

	str_del( id_string );
	str_del( selection_string );
	qp_qe_free( qe );
}


static void
do_nt_project ( void )
{
	char		*histo_string;
	char		*expr_string;
	char		*selection_string;
	char		*id_string, *id_path;
#ifdef LONG_LONG_CHAIN
	long long int	nevent, ifirst;
#else
	int		nevent, ifirst;
#endif
	int		i, err;
	int		id, idh, id_dim;
	SVec		expressions;
	QuerySrc *	qs;
	QueryExe *	qe;
	bool		labelx = FALSE, labely = FALSE;
	SMap		labelx_list, labely_list;

	labelx_list = 0;
	labely_list = 0;


	/* get all the parameters from kuip */
	histo_string = str_new( ku_getf() );
	expr_string = str_new( ku_getf() );
	selection_string = str_new( ku_getf() );
	C2FCBSTR( selection_string, PCCSE2.chcsel, 0 );
	nevent = ku_geti();
	ifirst = ku_geti();

	err = split_id_string( expr_string, &id_string, &expressions );
	if ( err != R_NOERR ) {
		str_del( expr_string );
		str_del( histo_string );
		str_del( selection_string );
		return;
	}

	/* load the histogram into memory */

	err = h_load_histo( histo_string, &idh, &id_dim );
	str_del( histo_string );

	if ( (err==R_NOERR) && (svec_entries(expressions) != id_dim) &&
		!h_flag_profile(idh) ) {
		sf_report( "Histogram %d is %d-dimensional, %d "
			"expression%s required instead of %d\n", idh, id_dim,
			id_dim, id_dim!=1?"s are":" is",
			svec_entries(expressions) );
		err = R_SYNTAX_ERROR;
	}
	if ( (err==R_NOERR) && h_flag_profile(idh) &&
		(svec_entries(expressions) != 2) ) {
		sf_report( "Histogram %d is a profile histogram,"
			"2 expressions required instead of %d\n",
			idh, svec_entries(expressions) );
		err = R_SYNTAX_ERROR;
	}

	/* load the ntuple into memory */

	if ( err == R_NOERR ) {
		err = h_load_nt( id_string, &id_path, &id );
	}
	str_del( id_string );

	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( expr_string );
		str_del( selection_string );
		svec_del( expressions );
		return;
	}

	qs = qp_qs_new( id_path, id, selection_string, expressions );
	str_del( id_path );
	svec_del( expressions );

	qe = qp_compile( qs, TRUE, &err );
	qp_qs_free( qs );

	if ( err != R_NOERR ) {
        	h_reset_dir();
		return;
	}

	if ( (id_dim == 2) && qe->expr_type[0] == D_STR ) {
		labely = TRUE;
		if ( ! HLABEQ( idh, "Y" ) ) {
			sf_report( "Histogram %d has no "
				"labels on Y axis\n", idh);
        		h_reset_dir();
			str_del( expr_string );
			str_del( selection_string );
			qp_qe_free( qe );
			return;
		}
		labely_list = h_get_labels( idh, "Y" );
	}

	if ( qe->expr_type[id_dim - 1] == D_STR ) {
		labelx = TRUE;
		if ( ! HLABEQ( idh, "X" ) ) {
			sf_report( "Histogram %d has no "
				"labels on X axis\n", idh);
        		h_reset_dir();
			str_del( expr_string );
			str_del( selection_string );
			qp_qe_free( qe );
			return;
		}
		labelx_list = h_get_labels( idh, "X" );
	}

	if ( HCPIAF.ntpiaf == CF_TRUE ) {
		int	istat;

		if ( ! NewPiaf() && qp_has_string_expr( qe ) ) {
			sf_report( "Character expressions not" 
				" available on old Piaf\n" );
        		h_reset_dir();
			str_del( expr_string );
			str_del( selection_string );
			qp_qe_free( qe );
			return;
		}
		qp_qe_free( qe );

		istat = 0;
		PFPING( 0, 1, istat );

		if ( istat == 0 ) {
			PFPUSH( istat );
		}

		if ( istat == 0 ) {
			PFHOUT( idh, istat );
		}

		if ( istat == 0 ) {
			char	*buf;
			int	len;

			/* rebuild command string */

			len = 45 + 1;
			len += strlen( expr_string );
			len += strlen( selection_string );

			buf = (char *) calloc( len, 1 ); qp_assert( buf != 0 );

#ifdef LONG_LONG_CHAIN
			sprintf( buf, "nt/proj%12d %s %s%12lld%12lld",
#else
			sprintf( buf, "nt/proj%12d %s %s%12d%12d",
#endif
				idh, expr_string, selection_string,
				nevent, ifirst );

			PFKUIP( buf, istat );

			free( (void *) buf );

		}

		QUEST.iquest[0] = istat;

	} else {
		QPCmd *		cmd;

		calc_event_range( qe->id, &ifirst, &nevent );

		cmd = qpcmd_new( id_dim == 1 && !h_flag_profile(idh) ?
				CMD_HFILL1 : CMD_HFILL2 );

		if ( h_flag_1d(idh) && !h_flag_profile(idh) ) {
			QPCmdHFill1	*h;

			/* 1 dim histograms */
			h = &cmd->u.hfill1;
			h->idh = idh;
			h->cvt_x = datatype_to_cvtcallback( qe->expr_type[0],
								labelx_list );

		} else {
			QPCmdHFill2	*h;

			/* 2 dim and profile histograms */
			h = &cmd->u.hfill2;
			h->idh = idh;
			h->cvt_x = datatype_to_cvtcallback( qe->expr_type[1],
								labelx_list );
			h->cvt_y = datatype_to_cvtcallback( qe->expr_type[0],
								labely_list );
		}

#ifdef PIAF_PROGRESS_METER
		PCNTWK.npass = 1;
		PCNTWK.ipass = 1;
		PCNTWK.misbyt = 0;
#endif
		qp_execute( qe, ifirst, nevent, cmd, &err );

		qpcmd_free( cmd );
		qp_qe_free( qe );
	}

	if ( labelx_list != 0 ) {
		smap_del( labelx_list );
	}
	if ( labely_list != 0 ) {
		smap_del( labely_list );
	}

	h_reset_dir();

	str_del( expr_string );
	str_del( selection_string );
}


static void
do_nt_plot ( void )
{
	char		*expr_string;
	char		*selection_string;
#ifdef LONG_LONG_CHAIN
	long long int	nevent, ifirst, nupd;
#else
	int		nevent, ifirst, nupd;
#endif
	char		*opt_string;
	int		idh;
	char		*id_string, *id_path;
	SVec		expressions;
	char		*hplot_option;
	int		i, err, id;
	bool		scatter_plot;
	QuerySrc *	qs;
	QueryExe *	qe;
	PlotOptions	opt, opt2;
	int		logx, logy, logz;

	/* get all the parameters from kuip */
	expr_string = str_new( ku_getf() );
	selection_string = str_new( ku_getf() );
	C2FCBSTR( selection_string, PCCSE2.chcsel, 0 );
	nevent = ku_geti();
	ifirst = ku_geti();
	nupd = ku_geti();
	opt_string = str_new( ku_getc() );
	idh = ku_geti();

	if ( HCPIAF.slavpf ) {
	  /* Old versions of PAW used to pass 0 as the default selection
	   * function. In this case 1. is substituted
	   */
	  if ( strcmp ( selection_string, "0" ) == 0) {
		str_del( selection_string );
		selection_string = str_new( "1." );
	  }
	}

	qp_plot_opt_scan( &opt, opt_string );

	/* work around coming from the "old" code (?) */
	if ( nevent == 0 ) {
		qp_abort( "No events selected ");
	}

	/* some conventions */

	if ( nevent < 0 ) {
		idh = -nevent;
		nevent = 99999999;
	}
	if ( ifirst < 0 ) {
		idh = -ifirst;
		ifirst = 1;
	}
	if ( nupd < 0 ) {
		idh = -nupd;
		nupd = 100000000;
	} else if ( nupd == 0 ) {
		nupd = 100000000;
	}

	err = split_id_string( expr_string, &id_string, &expressions );
	if ( err != R_NOERR ) {
		str_del( selection_string );
		str_del( opt_string );
		str_del( expr_string );
		return;
	}

	/* load the ntuple into memory */

	err = h_load_nt( id_string, &id_path, &id );
	str_del( id_string );
	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( selection_string );
		svec_del( expressions );
		str_del( opt_string );
		str_del( expr_string );
		return ;
	}

	qs = qp_qs_new( id_path, id, selection_string, expressions );
	str_del( id_path );
	svec_del( expressions );

	qe = qp_compile( qs, TRUE, &err );
	qp_qs_free( qs );

	if ( err != R_NOERR ) {
	        h_reset_dir();
		str_del( opt_string );
		str_del( selection_string );
		str_del( expr_string );
		return;
	}

	calc_event_range( qe->id, &ifirst, &nevent );

	/* filter out ntuple options  */
	opt2 = opt;
	opt2.S = FALSE;
	hplot_option = qp_plot_opt_gen( &opt2, TRUE );
	scatter_plot = ((*hplot_option == '\0')  || opt.line) && (!opt.N);
	str_del( hplot_option );

	logx = 0; logy = 0; logz = 0; /* bool might be converted, needs init */
	PAHLOG(logx,logy,logz);
	opt.logx = F2CLOGICAL(logx);
	opt.logy = F2CLOGICAL(logy);
	opt.logz = F2CLOGICAL(logz);

	/* reset the use flags in the buffer cache, free temp buffers */
	HBINIT1( 1 );

	if ( HCPIAF.ntpiaf == CF_TRUE ) {
		int	istat;

		if ( ! NewPiaf() && qp_has_string_expr( qe ) ) {
			sf_report( "Character expressions not"
				" available on old Piaf\n" );
			h_reset_dir();
			str_del( expr_string );
			str_del( selection_string );
			str_del( opt_string );
			qp_qe_free( qe );
			return;
		}

		istat = 0;
		PFPING( 0, 1, istat );

		if ( istat == 0 ) {
			PFPUSH( istat );
		}

		if ( istat == 0 ) {
			if ( HEXIST(idh) ) {
				if ( idh == 1000000 ) {
					HDELET( idh );
				} else {
					HRESET( idh, "" );
					PFHOUT( idh, istat );
				}
			}
		}

		if ( istat == 0 ) {
			char	*buf;
			int	len;

			/* rebuild command string */

			len = 128 + 1;
			len += strlen( expr_string );
			len += strlen( selection_string );
			len += strlen( opt_string );

			buf = (char *) calloc( len, 1 ); qp_assert( buf != 0 );

			/* NTUPLE/PLOT IDN [ UWFUNC NEVENT IFIRST NUPD OPTION IDH ] */
#ifdef LONG_LONG_CHAIN
			sprintf( buf, "nt/plot %s %s%12lld%12lld%12lld %s %12d",
#else
			sprintf( buf, "nt/plot %s %s%12d%12d%12d %s %12d",
#endif
				expr_string, selection_string,
				nevent, ifirst, nupd, 
				strlen(opt_string) > 0 ? opt_string : "!", idh );

			if ( qp_flags_get( "verbose" ) != 0 ) {
				sf_report( "PFKUIP(%s)\n",buf );
			}

			PFKUIP( buf, istat );

			free( (void *) buf );

		}

		QUEST.iquest[0] = istat;

	} else {

		switch( qe->nexpr ) {
		case 1:
			do_hplot_1d( qe, ifirst, nevent, idh, &opt, &err );
			break;

		case 2:
			if ( PROFILE_OPTION( opt ) ||
				(idh != 1000000 && h_flag_profile(idh) ) ) {
				do_prof_1d( qe, ifirst, nevent, idh,
					&opt, &err );
			} else if ( scatter_plot )  {
				do_scat_2d( qe, ifirst, nevent, idh,
					&opt, &err );
			} else {
				do_hplot_2d( qe, ifirst, nevent, idh,
					&opt, &err );
			}
			break;

		case 3:
			if ( opt.gouraud ) {
				do_gouraud( qe, ifirst, nevent, &opt, &err );

			} else {
				do_scat_3d( qe, ifirst, nevent, &opt, &err );
			}
			break;

		case 4:
			do_scat_4d( qe, ifirst, nevent, &opt, &err );
			break;

		default:
			sf_report( "Cannot plot %d-dimensional data\n",
				qe->nexpr );
			break;
		}

	}

	if ( err == R_NOERR && PCCSEL.ioptcn != 0 ) {
		PACSEL;
	}

	h_reset_dir();

	str_del( expr_string );
	str_del( selection_string );
	str_del( opt_string );
	qp_qe_free( qe );
}


static void
do_nt_gcut ( void )
{
	char		*cut_string;
	char		*expr_string;
	char		*selection_string;
#ifdef LONG_LONG_CHAIN
	int		cid;
	long long int	nevent, ifirst, nupd;
#else
	int		cid, nevent, ifirst, nupd;
#endif
	char		*opt_string;
	int		idh;
	char		*id_string, *id_path;
	SVec		expressions;
	char		*hplot_option;
	int		i, n, err, id, wkid;
	bool		scatter_plot;
	QuerySrc *	qs;
	QueryExe *	qe;
	PlotOptions	opt;
	int		logx, logy, logz;

	/* get all the parameters from kuip */

	cut_string = str_new( ku_getf() );
	cid = cut_get_cid( cut_string );
	str_del( cut_string );
	if ( cid == -1 ) {
		return;
	}

	expr_string = str_new( ku_getf() );
	selection_string = str_new( ku_getf() );
	C2FCBSTR( selection_string, PCCSE2.chcsel, 0 );
	nevent = ku_geti();
	ifirst = ku_geti();
	nupd = ku_geti();
	opt_string = str_new( ku_getc() );
	idh = ku_geti();
	wkid = ku_geti();

	qp_plot_opt_scan( &opt, opt_string );

	/* work around comming from the "old" code (?) */
	if ( nevent == 0 ) {
		qp_abort( "No events selected");
	}

	/* some conventions */

	if ( nevent < 0 ) {
		idh = -nevent;
		nevent = 99999999;
	}
	if ( ifirst < 0 ) {
		idh = -ifirst;
		ifirst = 1;
	}
	if ( nupd < 0 ) {
		idh = -nupd;
		nupd = 100000000;
	} else if ( nupd == 0 ) {
		nupd = 100000000;
	}

	err = split_id_string( expr_string, &id_string, &expressions );

	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( opt_string );
		str_del( expr_string );
		str_del( selection_string );
		return;
	}

	n = svec_entries( expressions );
	if ( n != 1 && n != 2 ) {
		sf_report( "A graphical cut can only be defined on one or two"
			" dimensional plots\n" );
		/* cleanup */
		str_del( opt_string );
		str_del( expr_string );
		str_del( selection_string );
		svec_del( expressions );
		return;
	}

	/* load the ntuple into memory */

	err = h_load_nt( id_string, &id_path, &id );
	str_del( id_string );
	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( opt_string );
		str_del( expr_string );
		str_del( selection_string );
		svec_del( expressions );
		return ;
	}


	qs = qp_qs_new( id_path, id, selection_string, expressions );
	str_del( id_path );
	svec_del( expressions );

	qe = qp_compile( qs, TRUE, &err );
	qp_qs_free( qs );

	if ( err != R_NOERR ) {
		h_reset_dir();
		str_del( opt_string );
		str_del( expr_string );
		str_del( selection_string );
		return;
	}

	calc_event_range( qe->id, &ifirst, &nevent );

	/* filter out ntuple options  */
	hplot_option = qp_plot_opt_gen( &opt, TRUE );
        if (opt.star  || opt.box   || opt.chr     || opt.col   || opt.cont  ||
            opt.lego  || opt.lego1 || opt.lego2   || opt.prof  || opt.profi ||
            opt.profs || opt.surf  || opt.surf1   || opt.surf2 || opt.surf3 ||
            opt.surf4 || opt.text  || opt.gouraud) {
	   scatter_plot = 0;
        } else {
	   scatter_plot = 1;
        }
	str_del( hplot_option );

	PAHLOG(logx,logy,logz);
	opt.logx = F2CLOGICAL(logx);
	opt.logy = F2CLOGICAL(logy);
	opt.logz = F2CLOGICAL(logz);

	/* reset the use flags in the buffer cache, free temp buffers */
	HBINIT1( 1 );

	if ( HCPIAF.ntpiaf == CF_TRUE ) {
		int	istat;

		if ( ! NewPiaf() && qp_has_string_expr( qe ) ) {
			sf_report( "Character expressions not"
				" available on old Piaf\n" );
			h_reset_dir();
			str_del( expr_string );
			str_del( selection_string );
			str_del( opt_string );
			qp_qe_free( qe );
			return;
		}		

		istat = 0;
		PFPING( 0, 1, istat );

		if ( istat == 0 ) {
			PFPUSH( istat );
		}

		if ( istat == 0 ) {
			if ( HEXIST(idh) ) {
				if ( idh == 1000000 ) {
					HDELET( idh );
				} else {
					PFHOUT( idh, istat );
				}
			}
		}

		if ( istat == 0 ) {
			char	*buf;
			int	len;

			/* rebuild command string */

			len = 128 + 1;
			len += strlen( expr_string );
			len += strlen( selection_string );
			len += strlen( opt_string );

			buf = (char *) calloc( len, 1 ); qp_assert( buf != 0 );

			/* NTUPLE/PLOT IDN [ UWFUNC NEVENT IFIRST NUPD OPTION IDH ] */
#ifdef LONG_LONG_CHAIN
			sprintf( buf, "nt/plot %s %s%12lld%12lld%12lld %s %12d",
#else
			sprintf( buf, "nt/plot %s %s%12d%12d%12d %s %12d",
#endif
				expr_string, selection_string,
				nevent, ifirst, nupd, 
				strlen(opt_string) > 0 ? opt_string : "!", idh );

			if ( qp_flags_get( "verbose" ) != 0 ) {
				sf_report( "PFKUIP(%s)\n",buf );
			}
			PFKUIP( buf, istat );

			free( (void *) buf );

		}

		QUEST.iquest[0] = istat;
	} else {

		switch( qe->nexpr ) {
		case 1:
			do_hplot_1d( qe, ifirst, nevent, idh, &opt, &err );
			break;

		case 2:
			if ( PROFILE_OPTION( opt ) ||
				(idh != 1000000 && h_flag_profile(idh) ) ) {
				do_prof_1d( qe, ifirst, nevent, idh,
					&opt, &err );
			} else if ( scatter_plot ) {
				do_scat_2d( qe, ifirst, nevent, idh,
					&opt, &err );
			} else {
				do_hplot_2d( qe, ifirst, nevent, idh,
					&opt, &err );
			}
			break;

		default:
			qp_abort( "Internal error\n" );
			break;
		}

	}

	if ( err == R_NOERR && PCCSEL.ioptcn != 0 ) {
		PACSEL;
	}

	str_del( selection_string );
	str_del( expr_string );
	str_del( opt_string );

	ku_alfa();
	if ( qe->nexpr == 2 && !PROFILE_OPTION( opt ) ) {
		float	xv[50], yv[50];
		int	n = 50;

		PAWLOC( n, xv, yv, -1, wkid, "-*" );
		if ( n < 3 ) {
			sf_report( "Need atleast 3 points for a two dimensional"
				" plot\n" );
		} else {
			gcut_add_2d( cid, qe->expr_str[0], qe->expr_str[1],
					n, xv, yv );
		}
	} else {
		float	xv[2], yv[2];
		int	n = 2;

		PAWLOC( n, xv, yv, -1, wkid, "-*" );
		if ( n < 2 ) {
			sf_report( "Need 2 points for a one dimensional"
				" plot\n" );
		} else {
			gcut_add_1d( cid, qe->expr_str[0], xv[0], xv[1] );
		}
	}

	h_reset_dir();

	qp_qe_free( qe );
}


static void
do_nt_dump ( void )
{
	char		*expr_string;
	char		*selection_string;
#ifdef LONG_LONG_CHAIN
	long long int	nevent, ifirst;
#else
	int		nevent, ifirst;
#endif
	char		*sep1_string;
	char		*sep2_string;
	char		*file_name;
	char		*id_string, *id_path;
	SVec		expressions;
	int		i, err, id;
	QuerySrc *	qs;
	QueryExe *	qe;
	FILE *		fp;

	/* get all the parameters from kuip */
	expr_string = str_new( ku_getf() );
	selection_string = str_new( ku_getf() );
	C2FCBSTR( selection_string, PCCSE2.chcsel, 0 );
	nevent = ku_geti();
	ifirst = ku_geti();
	file_name = str_new( ku_getf() );
	sep1_string = str_new( ku_getc() );
	sep2_string = str_new( ku_getc() );

	err = split_id_string( expr_string, &id_string, &expressions );
	str_del( expr_string );

	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( selection_string );
		str_del( sep1_string );
		str_del( sep2_string );
		str_del( file_name );
		return;
	}

	/* load the ntuple into memory */

	err = h_load_nt( id_string, &id_path, &id );
	str_del( id_string );
	if ( err != R_NOERR ) {
		/* cleanup */
		str_del( selection_string );
		str_del( sep1_string );
		str_del( sep2_string );
		str_del( file_name );
		svec_del( expressions );
		return ;
	}


	qs = qp_qs_new( id_path, id, selection_string, expressions );
	str_del( id_path );
	str_del( selection_string );
	svec_del( expressions );

	qe = qp_compile( qs, TRUE, &err );
	qp_qs_free( qs );

	if ( err != R_NOERR ) {
		/* cleanup */
		h_reset_dir();
		str_del( sep1_string );
		str_del( sep2_string );
		str_del( file_name );
		return;
	}

	calc_event_range( qe->id, &ifirst, &nevent );

	/* reset the use flags in the buffer cache, free temp buffers */
	HBINIT1( 1 );

	if ( strlen( file_name ) > 0 ) {
		extern int	errno;

		fp = fopen( file_name, "w" );
		if ( fp == 0 ) {
			sf_report( "Cannot open file %s (%s)\n", file_name,
				strerror( errno ) );
			h_reset_dir();
			str_del( sep1_string );
			str_del( sep2_string );
			str_del( file_name );
			return;
		}

	} else {
		fp = stdout;
	}

	if ( HCPIAF.ntpiaf == CF_TRUE ) {
		sf_report( "Ntuple on piaf not supported\n" );
	} else {
		QPCmd *		cmd;

		cmd = qpcmd_new( CMD_DUMP );
		cmd->u.dump.fp = fp;
		if ( strlen(sep1_string) != 0 ) {
			cmd->u.dump.sep1 = str_new(sep1_string);
		} else {
			cmd->u.dump.sep1 = str_new(" ");
		}
		if ( strlen(sep2_string) != 0 ) {
			cmd->u.dump.sep2 = str_new(sep2_string);
		} else {
			cmd->u.dump.sep2 = str_new(" ");
		}

#ifdef PIAF_PROGRESS_METER
		PCNTWK.npass = 1;
		PCNTWK.ipass = 1;
		PCNTWK.misbyt = 0;
#endif
		qp_execute( qe, ifirst, nevent, cmd, &err );

		if ( strlen( file_name ) > 0 ) {
			fclose( fp );
		}

		qpcmd_free( cmd );
	}

	h_reset_dir();

	qp_qe_free( qe );

	/* cleanup */
	str_del( sep1_string );
	str_del( sep2_string );
	str_del( file_name );

}


void
npantup_C( void )
{
	char		*cmd_path;
	char		*cmd_string;
	void		qp_pull_c_decl_obj( void );

	qp_pull_c_decl_obj();	/* get the common blocks declared */

	if ( setjmp( qp_abort_env ) != 0 ) {
		return;	/* we had a serious problem */
	} else {
		qp_abort_env_valid = 1;
	}

	cmd_path = str_new( ku_path() );
	cmd_string = strrchr( cmd_path, '/' );
	if ( cmd_string != 0 ) {
		cmd_string += 1;
	} else {
		cmd_string = cmd_path;
	}

	if ( strcasecmp( "loop", cmd_string ) == 0 ) {
		do_nt_loop();
	} else if ( strcasecmp( "project", cmd_string ) == 0 ) {
		do_nt_project();
	} else if ( strcasecmp( "plot", cmd_string ) == 0 ) {
		do_nt_plot();
	} else if ( strcasecmp( "dump", cmd_string ) == 0 ) {
		do_nt_dump();
	} else if ( strcasecmp( "gcut", cmd_string ) == 0 ) {
		do_nt_gcut();
	} else if ( strcasecmp( "scan", cmd_string ) == 0 ) {
		do_nt_scan();
	} else {
		sf_report( "*** Internal error: %s unknown\n", cmd_string );
	}

	str_del( cmd_path );

	/* just to be sure ... we do not want to come back here */
	qp_abort_env_valid = 0;
}


FCALLSCSUB0(npantup_C,PANNTU,panntu)	/* create fortran entry point */

