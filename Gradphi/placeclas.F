c------------------------------------------------------------------
c Places CLAS target in EVAC
c------------------------------------------------------------------
      subroutine placeclas
      implicit none
#     include "geometry.inc"
      integer nr,irot
      character*4 chmoth,chname,chonly
      real x,y,z
      chmoth='EVAC'
      chname='CLAS'
      nr=1
      x=0.
      y=0.
      z=cDist+cLong/2.
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
