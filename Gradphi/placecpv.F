c Places CPV in the HALL
      subroutine placecpv
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
      integer CPVR
c first we get rotational matrix
      CPVR=5000
      call ANGLES(cpvtheta,cpvphi,cpvpsi)
      call GSROTM(CPVR,thetax,phix,thetay,phiy,thetaz,phiz)            
      chname='CPV '
      chmoth='HALL'
      nr=1
      x=cpvx
      y=cpvy
      z=cpvz
      irot=CPVR
      chonly='ONLY'
      CALL GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
