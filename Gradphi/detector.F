c -----------------------------------------------------------------------------
c The volume containing the detector
c------------------------------------------------------------------------------      
      subroutine detector
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c box containing the detector
      chname='HALL'
      chshap='BOX '
      nmed=med_AIR
      par(1)=10000.
      par(2)=10000.
      par(3)=10000.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end
