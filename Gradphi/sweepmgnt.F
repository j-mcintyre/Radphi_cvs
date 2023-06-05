c----------------------------------------------------------------
c                sweeping magnet
c----------------------------------------------------------------
      subroutine sweepmgnt
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
      real x,y,z
c      
      chname='SMAG'
      chshap='BOX '
      nmed=med_magnetpole
      par(1)=mWide/2.
      par(2)=mHigh/2.
      par(3)=mLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='MGAP'
      chshap='BOX '
      nmed=med_magnetgap
      par(1)=mGap/2.
      par(2)=mHigh/2.
      par(3)=mLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='MGAP'
      chmoth='SMAG'
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end

