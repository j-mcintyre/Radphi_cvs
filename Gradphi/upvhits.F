*CMZ :          30/08/94  14.16.48  by  S.Ravndal
*-- Author : Andriy Kurylov
      subroutine upvhits
      implicit none
C.
c this subroutine saves hits for upv
C.
#define CERNLIB_TYPE 1
#include "geant321/gctmed.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcsets.inc"
#     include "geometry.inc"
#     include "hitconstants.inc"
C.
c reminder: numbers of UPVS blades
c		7
c
c	3		6
c	2		5
c	1		4
c		8
c
c looking upstream
c		
c we get here if only we are in the UPVS
c 
      real hi(7)
      real XD(3)
      save hi
      real dist
      integer IHIT
      integer iflag
C
      if (INWVOL.eq.1) then
c if entered a new volume (in this case a new copy of
c UPVS) or new track have to initialize hi.
        hi(1) = VECT(1)
        hi(2) = VECT(2)
        hi(3) = VECT(3)
        hi(4) = GETOT
        hi(5) = 0.
        hi(6) = TOFG !time of the entrance
        hi(7) = 0.
      else if (CHARGE.gt.0) then
c we want the actual tome of the hit and dist
        iflag=1 !XD is the position of the current point
                !inside the current daughter in LOCAL CS
        call GMTOD(VECT,XD,iflag) ! convert VECT to local frame
c dist-distance from the point of the hit to light-guide
        dist=lengthScint_UPV/2.0-XD(3)
        hi(7)=TOFG+dist/cEffectiveCPV ! delay time 
        hi(5) = hi(5)+DESTEP*exp(-cpvattenuation*dist)
      endif 
      if (INWVOL.eq.2.or.ISTOP.gt.0) then
        if (hi(5).gt.0) then
          call GSAHIT(ISET,IDET,ITRA,NUMBV,hi,IHIT)
        endif
      endif
      end
