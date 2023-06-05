c----------------------------------------------------------------
c                primary, secondary photon beam collimators
c----------------------------------------------------------------
      subroutine collim
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
c      
      chname='PCOL'
      chshap='TUBE'
      nmed=med_collimator
      par(1)=aBore1/2.
      par(2)=-1.	! as big as the mother
      par(3)=aLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='SCOL'
      chshap='TUBE'
      nmed=med_collimator
      par(1)=aBore2/2.
      par(2)=-1.	! as big as the mother
      par(3)=aLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      end

