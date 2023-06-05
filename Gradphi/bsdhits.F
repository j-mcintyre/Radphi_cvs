c This routine saves hits into JHITS from the BSD detector
c Written 6/7/99 by Tom Bogue
c this rountine assumes that it is only called when in BSR1,BSR2,or BSR3
c
      subroutine bsdhits
      implicit none
#define CERNLIB_TYPE 1
c
c include file defining 
c INWVOL, a flag that tells if you've entered a new volume
c VECT, the postion/momentum of the current particle
c DESTEP, the energy lost in the current step
c TOFG, the time of flight
#include "geant321/gctrak.inc"
c include file for volume parameters, including
c NAMES(NLEVEL) the name of the current volume
#include "geant321/gcvolu.inc"
c include file defining:
c IDET, ISET, and NUMBV
#include "geant321/gcsets.inc"
c include file for ITRA
#include "geant321/gckine.inc"
C.
      CHARACTER*4 CNAMES(20)
      EQUIVALENCE (NAMES(1),CNAMES(1))
C.
c
c include file for various hit parameters of the BSD, i.e.
c cEffectiveBSD, the effective speed of light in a scintillator
c cSubzEffectiveBSD, the effective speed of light in the z direction for
c  light traveling along a twisting scintillator ring
#     include "hitconstants.inc"
c the include file describing the geometry of the BSD, including:
c ScintHeight_BSD, the half height of the scintillators
#     include "geometry.inc"
c the hit parameters, which should be saved from call to call
      real hit
      dimension hit(7)
      save hit
c hit = (x,y,z,Ein,Elost,TOF,TDC)
c (x,y,z) is in local coordinates.  This is unique amongst the detectors 
c the time required for the signal from this hit to reach the PMT
      real dt,dE
c the initial position of the hit in local coordinates
      real x(3)
c the distance in the z- direction for the signal to 
c  travel to reach the PMT
      real dist
      integer IHIT
c
      if (INWVOL .eq. 1) then
c       initialize the system
c       note: VECT, and now hit, is in global coordinates
        hit(4)=GEKIN
        hit(5)=0
        hit(6)=TOFG
c now we need to calculate dt
c first step is to find the position in local coordinates
        CALL GMTOD(VECT,x,1)
c while we have the local coodinates, let's store them in (x,y,z)
        hit(1)=x(1)
        hit(2)=x(2)
        hit(3)=x(3)
c now to find the distance to the upstream end of the BSD 
c  (where the PMTs are)
        dist=x(3)+ScintHeight_BSD
c now to calculate dt
        if (CNAMES(NLEVEL) .eq. 'BSR1') then
          dt=dist/cSubzEffectiveRing1BSD
        elseif (CNAMES(NLEVEL) .eq. 'BSR2') then
          dt=dist/cSubzEffectiveRing2BSD
        else
          dt=dist/cSubzEffectiveRing3BSD
        endif
c now to find the time the signal reached the PMT
        hit(7)=TOFG+dt
      else
c       continue summing up the energy lost in this step
        dist=x(3)+ScintHeight_BSD
        if (CNAMES(NLEVEL) .eq. 'BSR1') then
          dE=DESTEP*exp(-dist*bsdAttenuationRing1)
        elseif (CNAMES(NLEVEL) .eq. 'BSR2') then
          dE=DESTEP*exp(-dist*bsdAttenuationRing2)
        else
          dE=DESTEP*exp(-dist*bsdAttenuationRing3)
        endif
        hit(5)=hit(5)+dE
c       and check to see if this is the last track for this particle 
c       this step
        if ((INWVOL .eq. 2) .or. (ISTOP .ne. 0)) then 
c         and if there was any energy lost to the detector
          if (hit(5) .gt. 0) then 
c           save this hit into JHITS
            call GSAHIT(ISET,IDET,ITRA,NUMBV,hit,IHIT)
          endif
        endif
      endif
      end
