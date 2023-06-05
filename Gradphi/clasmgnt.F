c----------------------------------------------------------------
c                post-CLAS sweeping magnet
c----------------------------------------------------------------
      subroutine clasmgnt
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
      real x,y,z
c      
      chname='WEEP'
      chshap='BOX '
      nmed=med_magnetpole
      par(1)=xWide/2.
      par(2)=xHigh/2.
      par(3)=xLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='WGAP'
      chshap='BOX '
      nmed=med_leadingap
      par(1)=xGap/2.
      par(2)=xHigh/2.
      par(3)=xLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='WHOL'
      chshap='TUBE'
      nmed=med_magnetgap
      par(1)=0.
      par(2)=xDiam/2.
      par(3)=xLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='WHOL'
      chmoth='WGAP'
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='WGAP'
      chmoth='WEEP'
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end

