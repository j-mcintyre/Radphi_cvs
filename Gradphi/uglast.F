*CMZ :          25/08/94  15.07.22  by  S.Ravndal
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE UGLAST
      implicit none
C
C Termination routine to print histograms and statistics
C
#define CERNLIB_TYPE 1
#include "geant321/gcflag.inc"
C
      CALL GLAST
C
C Save histograms to file
C
c     CALL HRPUT(0,'radphi.hist',' ')
C
c     CALL HISTDO(0)
C
      END
