c------------------------------------------------------------------
c Places WALL in the detector
c------------------------------------------------------------------
      subroutine placewall
      implicit none
#     include "geometry.inc"
      integer nr,irot
      character*4 chmoth,chname,chonly
      real x,y,z,shift
      shift=sDist+sThick/2.-wallz
      chmoth='HALL'
      chname='WALL'
      nr=1
      x=0.
      y=0.
      z=wallz
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
