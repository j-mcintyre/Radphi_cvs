c------------------------------------------------------------------
c Places BRAD in the vacuum pipe
c------------------------------------------------------------------
      subroutine placebrad
      implicit none
#     include "geometry.inc"
      integer nr,irot
      character*4 chmoth,chname,chonly
      real x,y,z
      chmoth='EVAC'
      chname='BRAD'
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
