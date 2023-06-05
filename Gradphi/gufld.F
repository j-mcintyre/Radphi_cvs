*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUFLD(X,F)
      implicit none
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to compute the magnetic field F             *
C.    *       at space point X                                         *
C.    *                                                                *
C.    *    ==>Called by : GHELIX,GRKUTA                                *
C.    *                                                                *
C.    ******************************************************************
C.
      REAL X(3),F(3)

#define CERNLIB_TYPE 1
#include "geant321/gcvolu.inc"
#include "geometry.inc"

C.
C.    ------------------------------------------------------------------
C.
      F(1)=0.
      F(2)=0.
      F(3)=0.
C
      if ((names(nlevel).eq.'MGAP').or.
     +    (names(nlevel).eq.'SMAG')) then
        F(1)=mField
      elseif ((names(nlevel).eq.'WGAP').or.
     +        (names(nlevel).eq.'WHOL').or.
     +        (names(nlevel).eq.'WMAG')) then
        F(1)=xField
      endif
      END
