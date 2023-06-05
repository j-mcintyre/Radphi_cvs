c------------------------------------------------------------------
c Places GTUB in the detector
c------------------------------------------------------------------
      subroutine placegtub
      implicit none
#     include "geometry.inc"
      integer nr,irot
      character*4 chmoth,chname,chonly
      real x,y,z,shift
      shift=sDist+sThick/2.-wallz
      chmoth='HALL'
      chname='GTUB'
      nr=1
      x=0.
      y=0.
      z=eLong+gLong/2. + 1.-shift	! 1cm gap between e-pipe and g-tube
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
