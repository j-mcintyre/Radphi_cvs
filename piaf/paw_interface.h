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
 *  paw_interface.h  --
 *	Declare interface to PAW.
 *
 *  Original: 12-Jan-1996 15:58
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: paw_interface.h,v 1.1.1.1 2004/03/20 20:05:13 radphi Exp $
 *
 *  $Log: paw_interface.h,v $
 *  Revision 1.1.1.1  2004/03/20 20:05:13  radphi
 *  cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
 *
 *  Revision 1.15  1999/11/02 15:37:41  couet
 *  - reorganise for NT
 *
 *  Revision 1.14  1999/07/05 16:13:03  couet
 *  - prototype added
 *
 *  Revision 1.13  1999/07/05 15:43:32  couet
 *  - hbook_interface.h in now replaced by hbook.h in CVSCOSRC
 *
 *  Revision 1.12  1999/06/28 15:08:54  couet
 *  - use now cfortran.h in $CVSCOSRC
 *
 *  Revision 1.11  1996/08/30 10:04:27  lecointe
 *  Restored Gouraud Shading in NT/PLOT
 *
 *  Revision 1.10  1996/08/21 12:55:33  lecointe
 *  Restore the spider plot in ntuple/scan
 *
 *  Revision 1.9  1996/05/24 09:16:06  dinofm
 *  Bug fixed when a FORTRAN selection function operates on a PIAF residing
 *  ntuple. The source file is sent to PIAF and compiled if '.f77' extension
 *  has been used.
 *
 *  Revision 1.8  1996/05/15 13:11:37  maartenb
 *  - Fix the CSELECT command.
 *
 *  Revision 1.7  1996/05/06 13:34:19  dinofm
 *  Code modified to take care of empty histograms detection on slave(s).
 *
 *  Revision 1.6  1996/04/23 18:38:08  maartenb
 *  - Add RCS keywords
 *
 *
 */

#ifndef CERN_PAW_INTERFACE
#define CERN_PAW_INTERFACE

#include	<cfortran/cfortran.h>

PROTOCCALLSFSUB1(FTNPRN,ftnprn,STRING)
#define FTNPRN(CHMESS) CCALLSFSUB1(FTNPRN,ftnprn,STRING,CHMESS)

PROTOCCALLSFSUB2(GETTP,gettp,PFLOAT,PFLOAT)
#define GETTP(THETA,PHI) CCALLSFSUB2(GETTP,gettp,PFLOAT,PFLOAT,THETA,PHI)

PROTOCCALLSFSUB2(HFIND,hfind,INT,STRING)
#define HFIND(IDD,CHROUT) CCALLSFSUB2(HFIND,hfind,INT,STRING,IDD,CHROUT)

PROTOCCALLSFSUB0(PACSEL,pacsel)
#define PACSEL CCALLSFSUB0(PACSEL,pacsel)

PROTOCCALLSFSUB3(GETNBINS,getnbins,PINT,PINT,PINT)
#define	GETNBINS(NX,NY,NZ) CCALLSFSUB3(GETNBINS,getnbins,PINT,PINT,PINT,NX,NY,NZ)

PROTOCCALLSFSUB3(PAHLOG,pahlog,PLOGICAL,PLOGICAL,PLOGICAL)
#define	PAHLOG(LOGX,LOGY,LOGZ) CCALLSFSUB3(PAHLOG,pahlog,PLOGICAL,PLOGICAL,PLOGICAL,LOGX,LOGY,LOGZ)

PROTOCCALLSFSUB9(PAPLOT,paplot,INT,STRING,STRING,INT,INT,INT,INT,INT,INT)
#define PAPLOT(ID,CHOPT,CHCASE,NUM,ICRANG,ICX1,ICX2,ICY1,ICY2) \
        CCALLSFSUB9(PAPLOT,paplot,INT,STRING,STRING,INT,INT,INT,INT,INT,INT,\
        ID,CHOPT,CHCASE,NUM,ICRANG,ICX1,ICX2,ICY1,ICY2)

PROTOCCALLSFSUB11(PASPI,paspi,INT,INT,STRING,FLOAT,INTV,INT,FLOATV,FLOATV,FLOATV,FLOATV,INT)
#define PASPI(ICHEVT,NVARS,CNAMES,RZONE,IVART,IZONE,CURRENT,LOW,HIGH,AVG,SPIDER_TYPE) \
        CCALLSFSUB11(PASPI,paspi,INT,INT,STRING,FLOAT,INTV,INT,FLOATV,FLOATV,FLOATV,FLOATV,INT,\
        ICHEVT,NVARS,CNAMES,RZONE,IVART,IZONE,CURRENT,LOW,HIGH,AVG,SPIDER_TYPE)

PROTOCCALLSFSUB8(PADRISO,padriso,INT,INT,INT,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV)
#define PADRISO(NX,NY,NZ,X,Y,Z,VALUE,S) \
        CCALLSFSUB8(PADRISO,padriso,INT,INT,INT,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,\
        NX,NY,NZ,X,Y,Z,VALUE,S) 

PROTOCCALLSFSUB1(PAUTIT,pautit,STRING)
#define PAUTIT(CHTIT) CCALLSFSUB1(PAUTIT,pautit,STRING,CHTIT)

PROTOCCALLSFSUB0(PAWCS,pawcs)
#define PAWCS CCALLSFSUB0(PAWCS,pawcs)

PROTOCCALLSFSUB4(PAWILD,pawild,STRING,STRING,INT,PINT)
#define PAWILD(CFIRST,CSECON,ILFIRS,IRESUL) CCALLSFSUB4(PAWILD,pawild,STRING,STRING,INT,PINT,CFIRST,CSECON,ILFIRS,IRESUL)

PROTOCCALLSFSUB6(PAWLOC,pawloc,PINT,FLOATV,FLOATV,INT,INT,STRING)
#define	PAWLOC(NP,XP,YP,NTPRI,IWKID,CHOPT) CCALLSFSUB6(PAWLOC,pawloc,PINT,FLOATV,FLOATV,INT,INT,STRING,NP,XP,YP,NTPRI,IWKID,CHOPT)

#ifdef LONG_LONG_CHAIN
  long long int pchevt_(char *NAME, int *LEN, int *ID, long long int *NEVT, int *IOP);
  long long int PCHEVT(char* name, int len, int id, long long int nevt, int iop);
#else
PROTOCCALLSFFUN5(INT,PCHEVT,pchevt,STRING,INT,INT,INT,INT)
#define	PCHEVT(NAME,LEN,ID,NEVT,IOP) CCALLSFFUN5(PCHEVT,pchevt,STRING,INT,INT,INT,INT,NAME,LEN,ID,NEVT,IOP)
#endif

PROTOCCALLSFSUB2(PCHNCD,pchncd,STRING,PINT)
#define	PCHNCD(PATH,IERR) CCALLSFSUB2(PCHNCD,pchncd,STRING,PINT,PATH,IERR)

#ifdef LONG_LONG_CHAIN
  void pcnext_(int *IDN, long long int *NCHROW, int *NDIM, int *NROW, int *IEND);
  void PCNEXT(int idn, long long int nchrow, int ndim, int nrow, int iend);
#else
PROTOCCALLSFSUB5(PCNEXT,pcnext,INT,PINT,PINT,PINT,PINT)
#define	PCNEXT(IDN,NCHROW,NDIM,NROW,IEND) CCALLSFSUB5(PCNEXT,pcnext,INT,PINT,PINT,PINT,PINT,IDN,NCHROW,NDIM,NROW,IEND)
#endif

PROTOCCALLSFSUB2(PFHOUT,pfhout,INT,PINT)
#define	PFHOUT(IDH,ISTAT) CCALLSFSUB2(PFHOUT,pfhout,INT,PINT,IDH,ISTAT)

PROTOCCALLSFSUB2(PFKUIP,pfkuip,STRING,PINT)
#define	PFKUIP(CHCMD,ISTAT) CCALLSFSUB2(PFKUIP,pfkuip,STRING,PINT,CHCMD,ISTAT)

PROTOCCALLSFSUB1(PFSOCK,pfsock,INT)
#define	PFSOCK(ISLAV) CCALLSFSUB1(PFSOCK,pfsock,INT,ISLAV)

PROTOCCALLSFSUB3(PFPING,pfping,INT,INT,PINT)
#define	PFPING(ISLAV,IACT,ISTAT) CCALLSFSUB3(PFPING,pfping,INT,INT,PINT,ISLAV,IACT,ISTAT)

PROTOCCALLSFSUB1(PFPUSH,pfpush,PINT)
#define	PFPUSH(ISTAT) CCALLSFSUB1(PFPUSH,pfpush,PINT,ISTAT)

PROTOCCALLSFSUB2(PFMINMAX,pfminmax,PFLOAT,PFLOAT)
#define	PFMINMAX(RMIN,RMAX) CCALLSFSUB2(PFMINMAX,pfminmax,PFLOAT,PFLOAT,RMIN,RMAX)

PROTOCCALLSFSUB2(PFLABELS,pflabels,PSTRING,PINT)
#define	PFLABELS(CBUF,LENBUF) CCALLSFSUB2(PFLABELS,pflabels,PSTRING,PINT,CBUF,LENBUF)

PROTOCCALLSFSUB1(PFEMPTY,pfempty,STRING)
#define	PFEMPTY(CBUF) CCALLSFSUB1(PFEMPTY,pfempty,STRING,CBUF)

PROTOCCALLSFSUB3(PFCSEX,pfcsex,INT,STRING,PINT)
#define	PFCSEX(LUNIN,CHFILE,ISTAT) CCALLSFSUB3(PFCSEX,pfcsex,INT,STRING,PINT,LUNIN,CHFILE,ISTAT)

#endif	/*	CERN_PAW_INTERFACE	*/
