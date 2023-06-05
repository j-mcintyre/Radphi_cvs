      block data backgroundbd
      implicit none
#     include "background.inc"
      data ebeam/5.65/
      data bgrate/0./
      data beamorigin/0.,0.,-3996./
      data bggate/-50.e-9,75.e-9/
      end

      subroutine superimposeBg
C
C Superimposes background in the form of random electrons incident on
C the bremsstrahlung target.  These electrons are forced in gukine() to
C commit a single act of discrete bremsstrahlung, and then die.  They
C are generated with a random time distribution over the duration of the
C ADC gate to simulate the actual conditions of a bremsstrahlung beam.
C
      implicit none
#define CERNLIB_TYPE 1
#include "geant321/gctrak.inc"
#     include "background.inc"
C
      integer nvert,nt
      integer ielec
      parameter (ielec=3)
      real dt
      real plab(3)
      real unif01(100)
      integer i,j

      if (bgrate.eq.0) return
      plab(1)=0
      plab(2)=0
      plab(3)=ebeam
      TOFG=bggate(1)+beamorigin(3)/3e10
      do i=1,99999
        call grndm(unif01,100)
        do j=1,100
          dt=-log(unif01(j))/bgrate
          TOFG=TOFG+dt*1e-9
          if (TOFG.gt.bggate(2)+beamorigin(3)/3e10) goto 10
          call gsvert(beamorigin,0,0,0,0,nvert);
          call gskine(plab,ielec,nvert,0,0,nt)
        enddo
      enddo
   10 continue
      end
