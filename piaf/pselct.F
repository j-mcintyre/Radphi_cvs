*
* R.T. Jones  July 10, 2002
* - Fixed the loop over bad slaves to save the names of the slaves
*   being suppressed before the array gets shifted away in psbads.F
*   To suppress this change, comment out the following line.
#define FIX_SLAVE_OVERWRITE_BUG 1
*
* $Id: pselct.F,v 1.1.1.1 2004/03/20 20:05:12 radphi Exp $
*
* $Log: pselct.F,v $
* Revision 1.1.1.1  2004/03/20 20:05:12  radphi
* cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
*
* Revision 1.1.1.1  1996/03/01 11:39:27  mclareni
* Paw
*
*
#include "paw/pilot.h"
*CMZ :  2.03/07 19/08/93  16.45.20  by  Fons Rademakers
*-- Author :    Fons Rademakers   16/08/93
      SUBROUTINE PSELCT(NSLAV, MASKT, ISLAV)
************************************************************************
*                                                                      *
* Make a select on NSLAV slave servers.                                *
*                                                                      *
* This routine returns via ISLAV the next socket from which something  *
* can be read.                                                         *
*                                                                      *
* *** Time out not yet implemented ***                                 *
* If ISLAV<0 then ABS(ISLAV) slaves timed out. The timed out slaves    *
* are considered dead and are added to the bad slave list. In case     *
* this happens a command should be re-executed. Slaves timeout after   *
* NSLAV times the time it took for the first slave to answer.          *
*                                                                      *
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
      MASK  = MASKT
*
      IF (MASK .EQ. 0) CALL PFMASK(NSLAV, MASK, -1)
*
      CALL PSMSEL(NSLAV, MASK, SLINPF, ISLAV, IBAD, ISTAT)
*
*-- In case some slaves timed out put them in the list of bad slaves
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
         ISLAV = -ISTAT
      ENDIF
*
      END
