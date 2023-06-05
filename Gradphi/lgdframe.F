c----------------------------------------------------------------
c                LGD frame
c----------------------------------------------------------------
      subroutine lgdframe
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c      
      chname='FRSD'
      chshap='BOX '
      nmed=med_IRON
      par(1)=frsdx/2
      par(2)=frsdy/2
      par(3)=frsdz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='FRHR'
      chshap='BOX '
      nmed=med_IRON
      par(1)=frhrx/2
      par(2)=frhry/2
      par(3)=frhrz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='FRTR'
      chshap='BOX '
      nmed=med_IRON
      par(1)=frtrx/2
      par(2)=frtry/2
      par(3)=frtrz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end

