c subroutine to place the BGV within the global volume (HALL)
c
      subroutine placebgv
      implicit none
c include file w/ nessesary geometrical data
#     include "geometry.inc"
c
c
c now just place the BGV into the global volume
      call GSPOS('BGV ',1,'HALL',0,0,Z_BGV,0,'ONLY')
      end
