*
* -R.T. Jones  July 30, 2003
*  Added support for ntuple chains longer than 2^31 events.
*  To do that we had to change a number of variables and arguments
*  from int to long long int, and hoped not to break anything.  So
*  far we believe that this should work on piaf without any changes
*  to paw itself.  To disable, comment out the following line.
*
#define LONG_LONG_CHAIN 1
*
* $Id: pawchn.inc,v 1.1.1.1 2004/03/20 20:05:13 radphi Exp $
*
* $Log: pawchn.inc,v $
* Revision 1.1.1.1  2004/03/20 20:05:13  radphi
* cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
*
* Revision 1.2  1996/03/01 16:49:56  cernlib
* All includes  had picked up a completely wrong #if defined(CERNLIB_UNIX), ie.
*  the real commons/declarations were missing on non-unix.
*
* Revision 1.1.1.1  1996/03/01 11:38:48  mclareni
* Paw
*
*
#ifndef CERNLIB_PAW_PAWCHN_INC
#define CERNLIB_PAW_PAWCHN_INC
*
*
* pawchn.inc
*
*CMZ :  2.02/01 10/05/93  00.52.38  by  Fons Rademakers
*-- Author :    Fons Rademakers   13/04/93
      LOGICAL         CHAIN
      CHARACTER*128   CFILE
#ifdef LONG_LONG_CHAIN
      COMMON /PAWCHN/ CHAIN, PADDING, NCHEVT, ICHEVT
      INTEGER*8       NCHEVT, ICHEVT
      INTEGER         PADDING
#else
      COMMON /PAWCHN/ CHAIN, NCHEVT, ICHEVT
#endif
      COMMON /PAWCHC/ CFILE
*
      CHARACTER*80    CURCHN
      COMMON /PAWCH2/ CURCHN
*
      LOGICAL         CHNXMT
      COMMON /PAWCH3/ CHNXMT
*

#endif
