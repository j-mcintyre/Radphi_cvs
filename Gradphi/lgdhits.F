c February 13, 2003
c   Added support for the registration of Cerenkov photon hits from
c   explicit Cerenkov following, instead of using the lead-glass
c   response model which is the default.
c - Richard Jones
c
c November 18, 2002
c   Modified the Cerenkov response in the LGD to take into account
c   the dependence of the N0 parameter on the angle of the charged
c   track in the shower with respect to the block axis.  There is
c   evidence from shower shape that this dependence is important.
c   To restore the old angle-independent behavior, enable this line.
c #define OLD_CERENKOV_N0 1
c - Richard Jones
c
*CMZ :          30/08/94  14.16.48  by  S.Ravndal
*-- Author : Andriy Kurylov
      subroutine lgdhits
      implicit none
C.
c this subroutine saves hits for lgd
C.
#define CERNLIB_TYPE 1
#include "geant321/gctmed.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcsets.inc"
#include "geant321/gconsp.inc"
#     include "hitconstants.inc"
#     include "geometry.inc"
C.
      real hi(6)
      save hi
      real cosch,sinch2,speed,length
      real zlocal
      real dNdx
      save dNdx
      integer IHIT
      real theta0
      integer itheta0
C
      if (IDTYPE.eq.1) then       ! indicates photocathode
        if ((IPART.eq.50).and.    ! indicates Cerenkov photon
     +      (DESTEP.gt.0)) then   ! photoelectron generated
          hi(1) = VECT(1)
          hi(2) = VECT(2)
          hi(3) = VECT(3)
          hi(4) = DESTEP
          hi(5) = 0.
          hi(6) = 1.
          call GSCHIT(ISET,IDET,ITRA,NUMBV,hi,3,IHIT)
        endif
      elseif (IDTYPE.eq.0) then   ! indicates LGD block
        if (INWVOL.eq.1) then
c if entered a new volume (in this case a new copy of
c BLOK) or new track have to initialize hi.
          hi(1) = VECT(1)
          hi(2) = VECT(2)
          hi(3) = VECT(3)
          if (NSTEP.gt.0) then
            hi(4) = GEKIN
          else
            hi(4) = 0. 
          endif  
          hi(5) = 0.
          hi(6) = 0. 
        else
          hi(6) = hi(6) + dNdx*STEP
        endif
        if (INWVOL.eq.2) then
          hi(5) = GEKIN
        endif
        if (INWVOL.eq.2.or.ISTOP.gt.0) then
          if (hi(6).gt.0) then
            call GSCHIT(ISET,IDET,ITRA,NUMBV,hi,3,IHIT)
          endif
        endif 
c now calculate rate of photons emitted per unit length for next step
c -speed of the particle
        speed = VECT(7)/GETOT
        cosch=1/(refindex*speed) ! cosine of the light cone angle
        if (cosch.le.1) then  ! now goes Cerenkov light
          sinch2=1-cosch**2   ! sin**2       
          zlocal=VECT(3)-lgdz
#if defined OLD_CERENKOV_N0
          dNdx = N0*sinch2
     +           *exp(zlocal*lgdattenuation)
#else
          if (abs(VECT(6)).lt.1) then
            theta0 = acos(VECT(6)) ! angle of track wrt block axis
          elseif (VECT(6).gt.0) then
            theta0 = 0
          else
            theta0 = acos(-1.)
          endif
          itheta0 = int(theta0*RADDEG) +1
          dNdx = cereN0(itheta0)*sinch2*(CHARGE**2)*cerewind *
     +           (lgdtail*exp((zlocal-lgdzmiddle)/lgdtaillen)
     +           +lgdrise*exp((zlocal-lgdzmiddle)/lgdriselen))
#endif
        else
          dNdx = 0
        endif
      endif
      end

      block data
#     include "hitconstants.inc"
      data cereN0/
     +     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     +     0.0000,0.0000,0.0000,0.0370,0.0740,0.1110,0.2220,
     +     1.0360,1.8130,2.4050,2.8490,3.2190,3.5150,3.8110,
     +     4.0330,4.1810,4.4030,4.5140,4.6250,4.7360,4.8470,
     +     4.9580,5.0320,5.0690,5.1060,5.1800,5.2170,5.1800,
     +     5.2170,5.2170,5.1800,5.1800,5.2170,5.1800,5.1060,
     +     5.1060,4.9950,4.9210,4.8470,4.7730,4.6990,4.6250,
     +     4.5510,4.5140,4.4770,4.3660,4.3660,4.3660,4.3290,
     +     4.2920,4.2920,4.2920,4.2180,4.1810,4.1810,4.1440,
     +     4.0700,4.0330,3.8850,3.8110,3.7370,3.6630,3.5520,
     +     3.4410,3.4410,3.2560,3.1450,3.0710,2.9230,2.8120,
     +     2.6640,2.5160,2.3680,2.1460,2.0350,1.8500,1.6650,
     +     1.4430,1.1840,0.9250,0.6290,0.3700,0.0740,90*0.0/
      end
