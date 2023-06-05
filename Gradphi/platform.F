c----------------------------------------------------------------
c                Platform for LGD and CPV
c----------------------------------------------------------------
      subroutine platform
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c      
      chname='PLA1'
      chshap='BOX '
      nmed=med_IRON
      par(1)=pla1x/2
      par(2)=pla1y/2
      par(3)=pla1z/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='PLA2'
      chshap='BOX '
      nmed=med_IRON
      par(1)=pla2x/2
      par(2)=pla2y/2
      par(3)=pla2z/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end

