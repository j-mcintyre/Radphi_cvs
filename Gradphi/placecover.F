c--------------------------------------------------------------
c Place the herculite cover
c---------------------------------------------------------------
      subroutine placecover
      implicit none
#     include "geometry.inc"
      character*4 chname, chmoth, chonly
      integer nr,irot
      real x,y,z
c
c front
c     
      chname='HFRO' 
      z=phfroz1
      x=phfrox1
      y=phfroy1
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='HFRO' 
      z=phfroz2
      x=phfrox2
      y=phfroy2
      chmoth='HALL'
      nr=2
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c      
c side
c            
      chname='HSID' 
      z=phsidz1
      x=phsidx1
      y=phsidy1
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='HSID' 
      z=phsidz2
      x=phsidx2
      y=phsidy2
      chmoth='HALL'
      nr=2
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c
c top
c            
      chname='HTOP' 
      z=phtopz1
      x=phtopx1
      y=phtopy1
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='HTOP' 
      z=phtopz2
      x=phtopx2
      y=phtopy2
      chmoth='HALL'
      nr=2
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
