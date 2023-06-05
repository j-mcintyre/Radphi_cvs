*CMZ :          30/08/94  14.16.48  by  S.Ravndal
*-- Author : Andriy Kurylov
      subroutine cpvhits
c
c this subroutine saves hits for cpv
c
      implicit none
#define CERNLIB_TYPE 1
#include "geant321/gctmed.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcsets.inc"
C.
      CHARACTER*4 CNAMES(20)
      EQUIVALENCE (NAMES(1),CNAMES(1))
C.
#     include "hitconstants.inc"
#     include "geometry.inc"
C.
c reminder: numbers of the CPV1 blades and their numbers
c   	7			15
c	5			13
c	3			11
c	1			9
c	2			10
c	4			12
c	6			14
c	8			16
c looking upstream
c
c numbers of the CPV2 blades
c	1			3
c	2			4
c looking upsrteam
c
c numbers of the CPV3 blades
c	4			9
c	2			7
c	1			6
c	3 			8
c	5			10
c looking upsrteam
c
c we get here if only we are in the CPV1, CPV2 or CPV3
c 
      real hi(7)
      real XD(3)
      save hi
      real dist
      integer IHIT
C
      if (INWVOL.eq.1) then
c if entered a new volume (in this case a new copy of
c CPV1, CPV2 or CPV3) or new track have to initialize hi.
        hi(1) = VECT(1)
        hi(2) = VECT(2)
        hi(3) = VECT(3)
        hi(4) = GETOT
        hi(5) = 0.
        hi(6) = TOFG !time of the entrance
        hi(7) = 0.
      else if (CHARGE .ne. 0) then
c we want the actual time of flight (meaning delay); 
c it depends on where we are
c but regardless we need local coordinates, so
        CALL GMTOD(VECT,XD,1)
        if (CNAMES(NLEVEL) .eq. 'CPSA') then
          dist=-XD(3)+LengthScint_CPVA/2.0
        elseif (CNAMES(NLEVEL) .eq. 'CPSB') then
          dist=-XD(3)+LengthScint_CPVB/2.0
        elseif (CNAMES(NLEVEL) .eq. 'CPSC') then
          dist=-XD(3)+LengthScint_CPVC/2.0
        elseif (CNAMES(NLEVEL) .eq. 'CPSD') then
          dist=-XD(3)+LengthScint_CPVD/2.0
        elseif (CNAMES(NLEVEL) .eq. 'CPSE') then
          dist=-XD(3)+LengthScint_CPVE/2.0
        else 
          print *,'error: unknown volume in cpvhits'
        endif
c now to calculate the time of signal arrival
        hi(7)=TOFG+dist/cEffectiveCPV
      hi(5) = hi(5)+DESTEP*exp(-cpvattenuation*dist)
      endif 
      if (INWVOL.eq.2.or.ISTOP.gt.0) then
        if (hi(5).gt.0) then
          call GSAHIT(ISET,IDET,ITRA,NUMBV,hi,IHIT)
        endif
      endif
      end
