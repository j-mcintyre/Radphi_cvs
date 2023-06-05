c Places RPD in the HALL
      subroutine placerpd
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
      integer RPDR
c first we get rotational matrix
      call ANGLES(rpdtheta,rpdphi,rpdpsi)
      RPDR=4000
      call GSROTM(RPDR,thetax,phix,thetay,phiy,thetaz,phiz)            
      chname='RPD '
      chmoth='HALL'      
      nr=1
      irot=RPDR
      chonly='ONLY'
      x=rpdx
      y=rpdy       
      z=rpdz
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
