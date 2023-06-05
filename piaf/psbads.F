*
* R.T. Jones  July 10, 2002
* - Fixed the loop that compresses the bad slaves out of the socket
*   info arrays, instead of always removing the first node from the list.
*   To disable this change comment out the following line.
#define FIX_SLAVE_OVERWRITE_BUG 1
*
* $Id: psbads.F,v 1.1.1.1 2004/03/20 20:05:12 radphi Exp $
*
* $Log: psbads.F,v $
* Revision 1.1.1.1  2004/03/20 20:05:12  radphi
* cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
*
* Revision 1.1.1.1  1996/03/01 11:39:27  mclareni
* Paw
*
*
#include "paw/pilot.h"
*CMZ :  2.04/02 16/11/93  15.08.48  by  Fons Rademakers
*-- Author :    Fons Rademakers   13/08/93
      SUBROUTINE PSBADS(HOST, ICLOSE)
************************************************************************
*                                                                      *
* Put a slave in the list of bad (dead) slaves. Optionally send a      *
* close message to the host.                                           *
*                                                                      *
************************************************************************
*
#include "hbook/hcdire.inc"
#include "hbook/hcpiaf.inc"
*
      CHARACTER*(*) HOST
*
*-- Find the slave
*
      DO 10 I = 1, NSLVPF
         IF (HOST .EQ. SLHOPF(I)) THEN
*
*-- Optionally close the host
*
            IF (ICLOSE .EQ. 1) THEN
               CALL PFSOCK(I)
               CALL CZCLOS(ISTAT)
               CALL PFSOCK(0)
            ENDIF
*
*-- Remove slave from tables
*
            NSLVPF = NSLVPF - 1
#ifdef FIX_SLAVE_OVERWRITE_BUG
            DO 20 J = I, NSLVPF
#else
            DO 20 J = 1, NSLVPF
#endif
               SLHOPF(J) = SLHOPF(J+1)
               SLINPF(J) = SLINPF(J+1)
               SLUTPF(J) = SLUTPF(J+1)
20          CONTINUE
*
            NCHTOP = NCHTOP - 1
#ifdef FIX_SLAVE_OVERWRITE_BUG
            DO 30 J = I, NCHTOP
#else
            DO 30 J = 1, NCHTOP
#endif
               CHTOP(J) = CHTOP(J+1)
               ICHTOP(J) = ICHTOP(J+1)
               ICHLUN(J) = ICHLUN(J+1)
               HFNAME(J) = HFNAME(J+1)
30          CONTINUE
*
*-- Add slave to list of bad slaves (so it will not be opened again)
*
            NSLBAD = NSLBAD + 1
            SLHBAD(NSLBAD) = HOST
            PRINT *,' Error accessing slave ',HOST(1:LENOCC(HOST)),
     +              ', slave closed'
*
            NSLAVE  = MIN(NSLAVE, NSLVPF)
*
         ENDIF
*
10    CONTINUE
*
      END
