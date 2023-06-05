c-----------------------------------------------------------------------
c    Target in the air highway
c-----------------------------------------------------------------------
      subroutine trgt
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
c other paremeters
      real x,y,z
c      

      chname='TRGT' 
      chshap='TUBE'
      nmed=med_BERYLLIUM_TGT
      par(1)=0.
      par(2)=tDiam/2
      par(3)=tLong/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c
      chname='WAY ' 
      chshap='BOX '
      nmed=med_AIR
      par(1)=2.
      par(2)=2.
      par(3)=0.95*(35.+(lgdz-upvz)/2) !35. is the half-thickness of LGD  
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c now putting the target
c first-its z-position
      z=lgdz+35.-par(3)
      z=targz-z+tLong/2.0
      x=0.
      y=0.
      irot=0
      nr=1
      chonly='ONLY'
      chname='TRGT'
      chmoth='WAY '
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c placing herculite cover inside the way
      z=lgdz+35.-par(3)
      z=phfroz1-z
      x=0.
      y=0.
      irot=0
      nr=1
      chonly='ONLY'
      chname='HFRO'
      chmoth='WAY '
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)                                                                              
      end
