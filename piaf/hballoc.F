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
* $Id: hballoc.F,v 1.1.1.1 2004/03/20 20:05:12 radphi Exp $
*
* $Log: hballoc.F,v $
* Revision 1.1.1.1  2004/03/20 20:05:12  radphi
* cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
*
* Revision 1.1.1.1  1996/01/16 17:07:56  mclareni
* First import
*
*
#include "hbook/pilot.h"
*CMZ :  4.22/09 22/07/94  12.48.52  by  Fons Rademakers
*-- Author :    Fons Rademakers   05/10/92
      SUBROUTINE HBALLOC(IDN,CHDIR,VAR,BLOCK,ITYPE,ISIZE,IFIRST,NELEM,
     +                   IBASE,IOFF,NUSE)
*.==========>
*.
*.   Dynamically allocate ISIZE*NELEM bytes for the storage of column
*.   VAR1 in BLOCK1 in Ntuple IDN. The allocated space can be accessed
*.   using IOFF, which is the offset from the base address IBASE.
*.   The first element of allocated space is IBASE(IOFF+1). If IOFF=0
*.   then there was not enough memory to allocate the space. When
*.   NUSE=1 then the buffer has been allocated for the first time and
*.   needs to be filled. If NUSE>1 then the buffer has already been
*.   used before and does not need to be filled again.
*.
*..=========> ( A.A.Rademakers )
*
      CHARACTER*(*)  CHDIR, VAR, BLOCK
#ifdef LONG_LONG_CHAIN
      INTEGER        IDN, ITYPE, ISIZE, IBASE(1), IBUF(1)
      INTEGER        IOFF, NUSE
      INTEGER*8      IFIRST, NELEM
#else
      INTEGER        IDN, ITYPE, ISIZE, NELEM, IBASE(1), IBUF(1)
      INTEGER        IOFF, NUSE, IFIRST
#endif
*
      LC = LENOCC(CHDIR)
      LV = LENOCC(VAR)
      LB = LENOCC(BLOCK)
*
      CALL HBALLO1(IDN, CHDIR, LC, VAR, LV, BLOCK, LB, ITYPE, ISIZE,
     +             IFIRST, NELEM, IBUF, NUSE)
*
      IF (NUSE .EQ. 0) THEN
         IOFF = 0
      ELSE
         IOFF = IBUF(1) - LOCF(IBASE(1))
      ENDIF
*
      END
