c--------------------------------------------------------------
c Place lgd frame
c---------------------------------------------------------------
      subroutine placeframe
      implicit none
#     include "geometry.inc"
      character*4 chname, chmoth, chonly
      integer nr,irot
      real x,y,z
c
c four vertical rods
c     
      chname='FRSD' 
      z=pfrsdz1
      x=pfrsdx1
      y=pfrsdy1
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRSD' 
      z=pfrsdz2
      x=pfrsdx2
      y=pfrsdy2
      chmoth='HALL'
      nr=2
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRSD' 
      z=pfrsdz3
      x=pfrsdx3
      y=pfrsdy3
      chmoth='HALL'
      nr=3
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRSD' 
      z=pfrsdz4
      x=pfrsdx4
      y=pfrsdy4
      chmoth='HALL'
      nr=4
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c      
c four horisontal rods
c            
      chname='FRHR' 
      z=pfrhrz1
      x=pfrhrx1
      y=pfrhry1
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRHR' 
      z=pfrhrz2
      x=pfrhrx2
      y=pfrhry2
      chmoth='HALL'
      nr=2
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRHR' 
      z=pfrhrz3
      x=pfrhrx3
      y=pfrhry3
      chmoth='HALL'
      nr=3
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRHR' 
      z=pfrhrz4
      x=pfrhrx4
      y=pfrhry4
      chmoth='HALL'
      nr=4
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c
c four transverse rods
c            
      chname='FRTR' 
      z=pfrtrz1
      x=pfrtrx1
      y=pfrtry1
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRTR' 
      z=pfrtrz2
      x=pfrtrx2
      y=pfrtry2
      chmoth='HALL'
      nr=2
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRTR' 
      z=pfrtrz3
      x=pfrtrx3
      y=pfrtry3
      chmoth='HALL'
      nr=3
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='FRTR' 
      z=pfrtrz4
      x=pfrtrx4
      y=pfrtry4
      chmoth='HALL'
      nr=4
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)      
      end
