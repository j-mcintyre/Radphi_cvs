c----------------------------------------------------------------
c Doughnut
c---------------------------------------------------------------
      subroutine dnut
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c
      chname='DNUT'
      chshap='TUBE'
      nmed=med_LEAD
      par(1)=dDiam1/2
      par(2)=dDiam2/2.
      par(3)=dLength/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end
