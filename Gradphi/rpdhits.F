*CMZ :          30/08/94  14.16.48  by  S.Ravndal
*-- Author : Andriy Kurylov
      subroutine rpdhits
      implicit none
C.
c this subroutine saves hits for rpd
C.
#define CERNLIB_TYPE 1
#include "geant321/gctmed.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcsets.inc"
C.
c we get here if only we are in the DSEN or THSN
c volumes that are sensitive volumes of the RPD
c 
      real hi(7)
      save hi
      integer IHIT
C
      if (INWVOL.eq.1) then
c if entered a new volume (in this case a new copy of
c DSEN or THSN) or new track have to initialize hi.
        hi(1) = VECT(1)
        hi(2) = VECT(2)
        hi(3) = VECT(3)
        hi(4) = GETOT
        hi(5) = 0.
        hi(6) = TOFG !time of the entrance
        hi(7) = hi(6) ! this is kept for compatibility with the subroutines
                      ! processing hits for CPV and UPV. This also gives 
                      ! the possibility to introuduce internad delay later  
      else if (abs(CHARGE).gt..1) then
        hi(5) = hi(5)+DESTEP
      endif 
      if (INWVOL.eq.2.or.ISTOP.gt.0) then
        if (hi(5).gt.0) then
          call GSAHIT(ISET,IDET,ITRA,NUMBV,hi,IHIT)
        endif
      endif
      end
