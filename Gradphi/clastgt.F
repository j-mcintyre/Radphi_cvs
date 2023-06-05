c----------------------------------------------------------------
c                CLAS target: CLAS
c----------------------------------------------------------------
      subroutine clastgt
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c      
      chname='CLAS'
      chshap='TUBE'
      nmed=med_CLAStarget
      par(1)=0.
      par(2)=cDiam/2.
      par(3)=cLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end

