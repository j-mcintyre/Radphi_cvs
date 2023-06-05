c--------------------------------------------------------------
c Place platform for the LGD and CPV
c---------------------------------------------------------------
      subroutine placeplatform
      implicit none
#     include "geometry.inc"
      character*4 chname, chmoth, chonly
      integer nr,irot
      real x,y,z
c     
      chname='PLA1' 
      z=ppla1z
      x=ppla1x
      y=ppla1y
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='PLA2' 
      z=ppla2z
      x=ppla2x
      y=ppla2y
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
