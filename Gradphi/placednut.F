c----------------------------------------------------------------
c Place Doughnut
c---------------------------------------------------------------
      subroutine placednut
      implicit none
#     include "geometry.inc"
      character*4 chmoth,chonly,chname
      integer nr,irot
      real x,y,z,shift
      shift=sDist+sThick/2-wallz
c      
      chmoth='HALL'
      chname='DNUT'
      nr=1
      x=0.
      y=0.
      z=wallz-sThick/2-dLength/2
      irot=0
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chmoth='GGAS'
      chname='DNUT'		! This is really the same object
      nr=2			! but embedding it in both volumes
      x=0			! allows them to overlap legally.
      y=0.
      z=z+shift-(eLong+gLong/2+1.)
      irot=0
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end      
