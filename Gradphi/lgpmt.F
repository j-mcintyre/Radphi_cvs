c----------------------------------------------------
c PMTs of the LGD
c----------------------------------------------------
      subroutine lgpmt
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      integer i,j,k
      real par(12)
      integer nt(28)
      character*4 chmoth,chonly
      integer nr,irot
c other parameters
      real x,y,z
      real zsense
c z shift of all stuff in box
      real z_Shift
      parameter (z_Shift=10.0)
c
      chname='PMTW' !box of air to contain PMTs
      chshap='BOX '
      nmed=med_AIR
      par(1)=blockx*Nx/2 !totally 28 rows; 13 in the upper half
      par(2)=blocky*Ny/2
      par(3)=5.08
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c now we have to specify the number of tubes per row
      nt(1)=10
      nt(2)=12
      nt(3)=16
      nt(4)=18
      nt(5)=20
      nt(6)=22
      nt(7)=24
      nt(8)=24
      nt(9)=26
      nt(10)=28
      nt(11)=28
      nt(12)=28
      nt(13)=28
      nt(14)=26            
c the lower half is symmetric
      do i=1,14
       nt(14+i)=nt(15-i)
      enddo
c defining the PMT
      chname='PMTL' !solid tube of pyrex glass
      chshap='TUBE'
      nmed=med_PYREX
      par(1)=0. 
      par(2)=1.58752
      par(3)=5.08
      zsense=-par(3)
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='PMTC' !photosensitive area on front of tube
      chshap='TUBE'
      nmed=med_cathode
      par(1)=0. 
      par(2)=1.50
      par(3)=0.10
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='PMTL'      
      nr=1
      irot=0
      chonly='ONLY'
      x=0.    
      y=0.
      z=zsense+par(3)
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c now we have to position the tubes
      chname='PMTL'
      chmoth='PMTW'
      irot=0
      chonly='ONLY'      
      nr=1      
      do i=1,28
       x=blockx*28/2-blockx*(i-1)-blockx/2
       z=0.      
       if(i.lt.14.or.i.gt.15) then
        do j=1,nt(i)
         y=blocky*(nt(i)-1)/2-blocky*(j-1)
         call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
         nr=nr+1
        enddo
       elseif(i.eq.14.or.i.eq.15) then
        do j=1,13
         do k =1,2
          y=(3-2*k)*(blocky*(nt(i)+1)/2-blocky*(j-1))
          call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
          nr=nr+1
         enddo
        enddo
       endif          
      enddo
c now positioning PMTW into the LGD 
      chname='PMTW'
      chmoth='LGD '      
      nr=1
      irot=0
      chonly='ONLY'
      x=0.
      y=0.
      z=28.00-z_Shift	! a few mm of air gap is sufficient
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
