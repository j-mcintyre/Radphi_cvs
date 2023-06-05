c----------------------------------------------------------------
c                Herculite Cover, Thickness 0.4mm
c----------------------------------------------------------------
      subroutine herccover
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c      
      chname='HFRO'
      chshap='BOX '
      nmed=med_PYREX
      par(1)=hfrox/2
      par(2)=hfroy/2
      par(3)=hfroz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='HSID'
      chshap='BOX '
      nmed=med_PYREX
      par(1)=hsidx/2
      par(2)=hsidy/2
      par(3)=hsidz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='HTOP'
      chshap='BOX '
      nmed=med_PYREX
      par(1)=htopx/2
      par(2)=htopy/2
      par(3)=htopz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end

