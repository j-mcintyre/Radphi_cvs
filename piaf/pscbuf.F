*
* R.T. Jones  July 10, 2002
* - Fixed the loop over bad slaves to save the names of the slaves
*   being suppressed before the array gets shifted away in psbads.F
*   To suppress this change, comment out the following line.
#define FIX_SLAVE_OVERWRITE_BUG 1
*
* $Id: pscbuf.F,v 1.1.1.1 2004/03/20 20:05:12 radphi Exp $
*
* $Log: pscbuf.F,v $
* Revision 1.1.1.1  2004/03/20 20:05:12  radphi
* cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
*
* Revision 1.1.1.1  1996/03/01 11:39:27  mclareni
* Paw
*
*
#include "paw/pilot.h"
*CMZ :  2.03/07 14/08/93  14.47.18  by  Fons Rademakers
*-- Author :    Fons Rademakers   14/08/93
      SUBROUTINE PSCBUF(NTOT, NSLAV, ISTAT)
************************************************************************
*                                                                      *
* Multicast NTOT characters from CHBUF to NSLAV slave servers.         *
*                                                                      *
* Equivalent semantics as CZPUTC.                                      *
* The array IBAD contains indices into the array SLINPF of dead        *
* slaves. The number of dead slaves is given by ISTAT.                 *
*                                                                      *
************************************************************************
*
#include "hbook/hcpiaf.inc"
*
      INTEGER        IBAD(MXSLPF)
#ifdef FIX_SLAVE_OVERWRITE_BUG
      character*100  mebad(mxslpf)
#endif
*
      ISTAT = 0
      IF (NSLAV .EQ. 0) RETURN
*
      CALL PSMBUF(NTOT, NSLAV, SLINPF, IBAD, ISTAT)
*
*-- In case some slave did not reply put it in the list of bad slaves
*
      IF (ISTAT .NE. 0) THEN
#ifdef FIX_SLAVE_OVERWRITE_BUG
         do i=1,istat
           mebad(i) = slhopf(ibad(i))
         enddo
#endif
         DO 10 I = 1, ISTAT
#ifdef FIX_SLAVE_OVERWRITE_BUG
            CALL PSBADS(mebad(i), 0)
#else
            CALL PSBADS(SLHOPF(IBAD(I)), 0)
#endif
10       CONTINUE
      ENDIF
*
      END
