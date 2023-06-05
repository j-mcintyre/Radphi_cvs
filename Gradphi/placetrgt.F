c LGD  in the HALL
      subroutine placetrgt
      implicit none
#     include "geometry.inc"
      character*4 chname, chmoth, chonly
      integer nr,irot
      real x,y,z
c     
      chname='WAY ' 
      z=lgdz+35.-0.95*(35.+(lgdz-upvz)/2)
      x=0.
      y=0.
      chmoth='HALL'
      nr=1
      chonly='ONLY'
      irot=0
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
