c subroutine to place the BSD into the main volume ('HALL')
c
      subroutine placebsd
      implicit none
#     include "geometry.inc"
c
c
c Now to place the whole BSD into the global volume
      call GSPOS('BSD ',1,'HALL',0,0,Z_BSD,0,'ONLY')
      end
