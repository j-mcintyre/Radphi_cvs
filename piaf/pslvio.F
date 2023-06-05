*
* R.T. Jones  July 13, 2002
* - Suppressed the escape clauses (goto 20) upon communications errors with
*   slaves.  If there are any slaves left then you must catch their replies
*   or you will lose sync and the slaves will crash on the next command.
*   To suppress this fix, comment out the following line.
#define FIX_STATUS_SYNC_BUG 1
*
* $Id: pslvio.F,v 1.4 2005/11/04 15:00:07 radphi Exp $
*
* $Log: pslvio.F,v $
* Revision 1.4  2005/11/04 15:00:07  radphi
* pslvio.F
*    - trap internal write errors with write(str,fmt,err=NNN) that appear
*      whenever a slave reports statistics with overflow values "***" in one
*      of the fields [rtj]
* qp_execute.c
*    - report an error message when a read error occurs during ntuple
*      processing; prior behavior was to abort processing and return the
*      results up to the error without a peep what happened! [rtj]
*
* Revision 1.3  2005/11/01 16:31:11  radphi
* psexec.F, pslvio.F
*    - changed the field lengths on piaf/status messages to prevent
*      the 80-character message from overflowing [rtj]
*
* Revision 1.2  2005/11/01 15:34:30  radphi
* psexec.F
*    - just a spelling error in the comments [rtj]
* pslvio.F
*    - make the Format statement more robust [rtj]
* qp_execute.c, qp_hbook_if.c qp_hbook_if.h
*    - fix a bug in qp_exec_init (qp_execute.c) where it was passing a short
*      row count (PCNEXT_OPTIMIZATION only, in the case where the slave
*      computes a non-contiguous slice of the ntuple chain) and so reserving
*      too little space in the ntuple data cache to hold the data [rtj]
*    - introduced a new function h_nrows_piaf() to return the TOTAL number
*      of rows in a slave's ntuple processing slice, as opposed to the size
*      of the first contiguous segment which is returned by h_range_piaf()
*      and h_next_nt(). [rtj]
*
* Revision 1.1.1.1  2004/03/20 20:05:12  radphi
* cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
*
* Revision 1.1.1.1  1996/03/01 11:39:28  mclareni
* Paw
*
*
#include "paw/pilot.h"
*CMZ :  2.05/08 09/06/94  18.58.54  by  Fons Rademakers
*-- Author :    Fons Rademakers   18/08/93
      SUBROUTINE PSLVIO(NSLAV, ISTAT)
************************************************************************
*                                                                      *
* Get the real and virt. memory IO statistics from the slave servers.  *
*                                                                      *
************************************************************************
*
#include "hbook/hcpiaf.inc"
*
      SRIN   = 0
      SROUT  = 0
      SVMIO  = 0
      STUSER = 0
      STCPU  = 0
      ISTAT  = 0
      IF (NSLAV .EQ. 0) RETURN
*
*-- Ask I/O and CPU status from slaves
*
      CALL PSCAST('IOSTAT', NSLAV, ISTAT)
      IF (ISTAT .NE. 0) THEN
         CALL PIAFLG(2,'Error broadcasting IOSTAT')
#ifndef FIX_STATUS_SYNC_BUG
         GOTO 20
#endif
      ENDIF
*
*-- Loop over slave server replies in the order they are coming in
*
      CALL PFMASK(NSLAV, SMASK, -1)
10    CALL PSELCT(NSLAV, SMASK, ISLAV)
      IF (ISLAV .LT. 0) THEN
         ISTAT = 1
         CALL PIAFLG(2,'Slaves do not respond to IOSTAT')
#ifndef FIX_STATUS_SYNC_BUG
         GOTO 20
#endif
      ENDIF
*
*-- Get slaves real and virtual IO statistics
*
      CALL PFSOCK(ISLAV)
      CALL CZGETA(CHSMPF, ISTAT)
      IF (ISTAT .NE. 0) THEN
         CALL PIAFLG(2,'Error getting I/O statistics')
         GOTO 15
      ENDIF
*
      READ(CHSMPF,'(3F15.0,2F15.3)',ERR=12) RIN, ROUT, VIO, USER, CPU
      SRIN   = SRIN   + RIN
      SROUT  = SROUT  + ROUT
      SVMIO  = SVMIO  + VIO
      STUSER = STUSER + USER
      STCPU  = STCPU  + CPU
12    CONTINUE   ! no heartbreak if overflows appear in some fields
*
15    CALL PFMASK(ISLAV, SMASK, 0)
      IF (SMASK .EQ. 0) GOTO 20
      GOTO 10
*
20    CALL PFSOCK(0)
*
**    WRITE(CHSMPF, '(A,2F16.0,I2)') 'PSLVIO:', SRIN, SROUT, NSLAV
**    CALL PIAFLG(1,CHSMPF(1:LENOCC(CHSMPF)))
*
      END
