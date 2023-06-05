c----------------------------------------------------------------
c                Electron beam pipe
c----------------------------------------------------------------
      subroutine beampipe
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
      real x,y,z
c      
      chname='EPIP'
      chshap='TUBE'
      nmed=med_ebeampipe
      par(1)=0.
      par(2)=eDiam/2.+eThick
      par(3)=eLong
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='EVAC'
      chshap='TUBE'
      nmed=med_VACUUM
      par(2)=eDiam/2.
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='WIND'
      chshap='TUBE'
      nmed=med_vacwindow
      par(1)=0.
      par(2)=-1.
      par(3)=wThick/2.
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='EPIP'
      chname='EVAC'
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chmoth='EVAC'
      chname='WIND'
      nr=1
      z=eLong-wThick/2.
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end

