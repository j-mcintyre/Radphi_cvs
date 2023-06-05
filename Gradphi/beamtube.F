c----------------------------------------------------------------
c                Photon beam tube
c----------------------------------------------------------------
      subroutine beamtube
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
      real x,y,z
c      
      chname='GTUB'
      chshap='TUBE'
      nmed=med_gbeamtube
      par(1)=0.
      par(2)=gDiam/2.+gThick
      par(3)=gLong/2.
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='GGAS'
      chshap='TUBE'
      nmed=med_gtubegas
      par(2)=gDiam/2.
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='GTUB'
      chname='GGAS'
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chmoth='GGAS'
      chname='WIND'
      nr=2
      z=(gLong-wThick)/2.
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=3
      z=-z
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
