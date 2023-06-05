c ------------------------------------------------------------ 
c LGD
c--------------------------------------------------------------
      subroutine lgd
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      integer nxdiv,nydiv
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
c angles for rotational matrices      
      real thetax,thetay,thetaz,phix,phiy,phiz 
      common thetax,thetay,thetaz,phix,phiy,phiz      
c other paremeters
      real x,y,z
c shift of box
      real z_Shift
      parameter (z_Shift=10.0)
c-------------------------------------------------------------------------------             
c lead glass wall-general containing volume
      chname='LGD '
      chshap='BOX '
      nmed=med_AIR
      par(1)=blockx*Nx/2
      par(2)=blocky*Ny/2
      par(3)=33.2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='WAL1' !box to contain blocks left and right of the hole
      chshap='BOX '
      nmed=med_wrap
      par(1)=blockx*(Nx-2)/4 !totally 28 rows; 13 in each wing
      par(2)=blocky*Ny/2
      par(3)=blockz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c divide WAL1 into blocks
c these blocks are sensitive detectors
c only some of them have tubes but this will be accounted for 
c when saving hits; now they are equal
      nxdiv=(Nx-2)/2
      nydiv=Ny         
      call GSDVN('GBX1','WAL1',nxdiv,1)
      call GSDVN('LGB1','GBX1',nydiv,2)
c now have to create insertion volumes
      chname='WAL2' !box to contain blocks above/below the hole
      chshap='BOX '
      nmed=med_wrap
      par(1)=blockx          !totally 2 columns; 13 in each column
      par(2)=blocky*(Ny-2)/4
      par(3)=blockz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c lets divide them into smaller blocks
      nxdiv=2
      nydiv=(Ny-2)/2         
      call GSDVN('GBX2','WAL2',nxdiv,1)
      call GSDVN('LGB2','GBX2',nydiv,2)
c now we have to place actual blocks in the divisions
      chname='BLKC' !lead glass block with loose wrapping
      chshap='BOX ' 
      nmed=med_gap
      par(1)=blockx1/2 + 100e-4
      par(2)=blocky1/2 + 100e-4
      par(3)=blockz/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='BLOK' !actual lead glass block
      chshap='BOX ' 
      nmed=med_LEAD_GLASS_SF3
      par(1)=blockx1/2 
      par(2)=blocky1/2
      par(3)=blockz1/2
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)    
c now we place the BLOK's into cells
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0.
      chonly='ONLY'
      chname='BLOK'
      chmoth='BLKC' 
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=1
      x=0.
      y=0.
      z=0.
      irot=0.
      chonly='ONLY'
      chname='BLKC'
      chmoth='LGB1' 
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=2
      chmoth='LGB2'
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c now we place two copies of WAL1 and WAL2 into the unit LGD 
      chname='WAL1'
      nr=1
      x=blockx*(Nx+2)/4
      y=0.
      z=-z_Shift
      irot=0.
      chonly='ONLY'
      chmoth='LGD ' 
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c now second copy of WAL1     
      nr=2
      x=-x
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c now the first copy of WAL2
      chname='WAL2'
      nr=1
      x=0.
      y=blocky*(Ny+2)/4
      z=-z_Shift
      irot=0.
      chonly='ONLY'
      chmoth='LGD ' 
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c now second copy of WAL2    
      nr=2
      y=-y
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      end
