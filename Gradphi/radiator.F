c----------------------------------------------------------------
c                bremsstrahlung radiator: BRAD
c----------------------------------------------------------------
      subroutine radiator
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c      
      chname='BRAD'
      chshap='TUBE'
      nmed=med_radiator
      par(1)=0.
      par(2)=rDiam/2.
      par(3)=rThick/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end

