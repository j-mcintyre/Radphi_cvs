*CMZ :          30/08/94  16.03.40  by  S.Ravndal
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE UGINIT
      implicit none
C
C To initialise GEANT/USER  program and read data cards
C
#define CERNLIB_TYPE 1
#include "geant321/gckine.inc"
#include "geant321/gcflag.inc"
#     include "geometry.inc"
#     include "background.inc"
C
C Open user files
C
C     CALL UFILES
C
C Initialise GEANT
C
      CALL GINIT
C
C Define user FFREAD data cards (format free input)
C
      CALL FFKEY('bgrate',bgrate,1,'REAL')
C
C Read the data cards
C
      OPEN(UNIT=4,FILE='control.in',STATUS='UNKNOWN')
      CALL FFSET('LINP',4)
      CALL GFFGO
C
C Initialise Zebra structure
C
      CALL GZINIT
C
C Geometry and materials description
C
      CALL major
C
C Print the defined materials, tracking media and volumes
C
C      CALL GPRINT('MATE',0)
C      CALL GPRINT('TMED',0)
C      CALL GPRINT('VOLU',0)
C
C Define user histograms and n-tuple
C
      CALL UHINIT
C
      END
