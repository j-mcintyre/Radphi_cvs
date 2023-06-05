c Places UPV in the HALL
      subroutine placeupv
      implicit none
#     include "geometry.inc"
      real x,y,z
      x=upvx
      y=upvy
      z=upvz		   
      call GSPOS('UPV ',1,'HALL',x,y,z,0,'ONLY')
      end
