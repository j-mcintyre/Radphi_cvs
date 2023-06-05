c Places LGD in the HALL
      subroutine placelgd
      implicit none
#     include "geometry.inc"
      character*4 chname, chmoth, chonly
      integer nr,irot
c angles for rotational matrices      
      real thetax,thetay,thetaz,phix,phiy,phiz 
      common thetax,thetay,thetaz,phix,phiy,phiz      
c other paremeters
      real x,y,z
c these are rotation matrices
      integer LGDR
c first we get rotational matrix
      call ANGLES(lgdtheta,lgdphi,lgdpsi)
      LGDR=1
      call GSROTM(LGDR,thetax,phix,thetay,phiy,thetaz,phiz)            
      chname='LGD '
      chmoth='HALL'      
      nr=1
      irot=LGDR
      chonly='ONLY'
      x=lgdx  
      y=lgdy       
      z=lgdz
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
