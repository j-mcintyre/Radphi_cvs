c----------------------------------------------------------------
c                CLAS magnet radiation shield
c----------------------------------------------------------------
      subroutine shield
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c      
      chname='HOLE'
      chshap='CONE'
      nmed=med_CLASshield
      par(1)=oLong/2.
      par(2)=oDiam1/2.
      par(3)=oThick
      par(4)=oDiam2/2.
      par(5)=oThick
      npar=5
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end

