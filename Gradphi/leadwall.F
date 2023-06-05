c----------------------------------------------------------------
c                upstream RadPhi lead shielding wall
c----------------------------------------------------------------
      subroutine leadwall
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      integer nr,irot
      character*4 chmoth,chonly
      real x,y,z
      real par(12)
c      
      chname='WALL'
      chshap='BOX '
      nmed=med_leadbrick
      par(1)=sWide/2
      par(2)=sHigh/2
      par(3)=sThick/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='APER'
      chshap='TUBE'
      nmed=med_AIR
      par(1)=0.
      par(2)=sDiam/2
      par(3)=sThick/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='WALL'
      chname='APER'
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)            
      end

