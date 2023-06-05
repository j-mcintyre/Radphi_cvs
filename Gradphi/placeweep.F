c-----------------------------------------------------------------
c WEEP in the detector
c-----------------------------------------------------------------
      subroutine placeweep
      implicit none
#     include "geometry.inc"
      integer nr,irot
      character*4 chmoth,chname,chonly
      real x,y,z
      chmoth='GGAS'
      chname='WEEP'
      nr=3
      x=0.
      y=0.
      z=oLong+xLong/2.-gLong/2.+50.
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
