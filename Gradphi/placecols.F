c------------------------------------------------------------------
c Places collimators in the detector
c------------------------------------------------------------------
      subroutine placecols
      implicit none
#     include "geometry.inc"
      integer nr,irot
      character*4 chmoth,chname,chonly
      real x,y,z
      chmoth='EVAC'
      chname='PCOL'
      nr=1
      x=0.
      y=0.
      z=aDist1+aLong/2.
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chmoth='EVAC'
      chname='SCOL'
      nr=1
      x=0.
      y=0.
      z=aDist2+aLong/2.
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chmoth='EVAC'
      chname='SMAG'
      nr=1
      x=0.
      y=0.
      z=(aLong+aDist1+aDist2)/2. - mlong
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=2
      x=0.
      y=0.
      z=(aLong+aDist1+aDist2)/2. + mlong
      irot=0     
      chonly='ONLY'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
