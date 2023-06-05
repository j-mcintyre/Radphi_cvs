c-----------------------------------------------------------------------------------                      
c            RECOIL PROTON DETECTOR
c-------------------------------------------------------------------------------
      subroutine rpd
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
c angles for rotational matrices      
      real thetax,thetay,thetaz,phix,phiy,phiz 
      common thetax,thetay,thetaz,phix,phiy,phiz      
c other paremeters
      real x,y,z,x1,y1,delz
c these are rotation matrices      
      integer PMTDR, DSENR, LGTHR1, LGTHR2
      integer THSNR1, THSNR2, PMTER1, PMTER2
c------------------------------------------------------------------------------
      delz=rpdz*2
c      
c this detector consists of two pars each of which has symetry of a polygon
c let's use this symetry
c first goes the front par-hexagon
      chname='RPD ' 
      chshap='PGON'
      nmed=med_AIR
      par(1)=0.
      par(2)=360. !whole circle
      par(3)=6 !hexagon
      par(4)=2 !cross-section is constant everywhere
      par(5)=-66.9/2
      par(6)=6.
      par(7)=110./2
      par(8)=66.9/2
      par(9)=6.
      par(10)=110./2 
      npar=10
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c now we have to divide this thing in 6 divisions in phi direction
      call GSDVN('CEL1','RPD ',6,2)
c it is enough to place pars of RPD in one of the CELs
c now local x-axes passes through the center of the division
c define the pars: sensitive vpolume DSEN, scintillator LGDS and the tube 
      chname='LGDS'
      chshap='TRD2'
      nmed=med_ACRYLIC
      par(1)=2.54
      par(2)=0.3175
      par(3)=2.54
      par(4)=10.1340
      par(5)=7.62
      npar=5
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)      
      chname='DSEN'
      chshap='TRAP'
      nmed=med_PLASTIC_SCINT
      par(1)=0.318
      par(2)=0.0
      par(3)=0.0
      par(4)=3.537
      par(5)=10.014
      par(6)=7.367
      par(7)=0.0
      par(8)=3.537
      par(9)=10.254
      par(10)=7.607
      par(11)=0.0
      npar=11
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='PMTD'
      chshap='TUBE'
      nmed=med_PYREX
      par(1)=3.0
      par(2)=3.175
      par(3)=8.89
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c now we put all these volumes in the division CEL1
c with the right rotation matrices
c rotmat for DSEN      
      DSENR=1000
      thetax=90.  
      phix=90.
      thetay=40.
      phiy=180.
      thetaz=50.
      phiz=0.
      call GSROTM(DSENR,thetax,phix,thetay,phiy,thetaz,phiz)
c the following is the rottation matrix for LGDS relative to CEL1
      PMTDR=2000
      thetax=50.
      phix=0.
      thetay=90.
      phiy=90.
      thetaz=40.
      phiz=180.
      call GSROTM(PMTDR,thetax,phix,thetay,phiy,thetaz,phiz)      
c rotmat for LGDSR is the same as for PMTDR
c now positioning these elements into CEL1
      chname='DSEN'
      chmoth='CEL1'      
      nr=1
      irot=DSENR
      chonly='ONLY'
#ifdef SOLARIS
      x=1.275*15.639*SIND(50.)     ! these are hardwired in here 
      y=0.0			   ! just for the time being  
      z=1.275*15.639*COSD(50.)-delz/2
#else
      x=1.275*15.639*DSIN(50.d0)   ! these are hardwired in here 
      y=0.0			   ! just for the time being  
      z=1.275*15.639*DCOS(50.d0)-delz/2
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='LGDS'
      chmoth='CEL1'      
      nr=1
      irot=PMTDR
      chonly='ONLY'
#ifdef SOLARIS
      x=22.8699*SIND(79.32)        ! see the comment above
      y=0.0
      z=22.8699*COSD(79.32)-delz/2      
#else
      x=22.8699*DSIN(79.32d0)      ! see the comment above
      y=0.0
      z=22.8699*DCOS(79.32d0)-delz/2      
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='PMTD'
      chmoth='CEL1'      
      nr=1
      irot=PMTDR
      chonly='ONLY'
#ifdef SOLARIS
      x=34.1330*SIND(104.25)	   ! see the comment for DSEN
      y=0.0
      z=34.1330*COSD(104.25)-delz/2
#else
      x=34.1330*DSIN(104.25d0)	   ! see the comment for DSEN
      y=0.0
      z=34.1330*DCOS(104.25d0)-delz/2
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
c it is enough to place pars of RPD-rear par in ther CEL1
c there will be two pars of the rear par in one cell-because the rear par 
c contains 12 scintillators 
c now local x-axes passes through the center of the division
c define the pars: sensitive vpolume THSN, scintillator LGTH and the tube PTME 
      chname='LGTH'
      chshap='TRAP'
      nmed=med_ACRYLIC
      par(1)=9.525
      par(2)=0.0
      par(3)=0.0
      par(4)=5.795
      par(5)=7.
      par(6)=5.005
      par(7)=0.
      par(8)=2.54
      par(9)=2.54
      par(10)=1.66557 !fixed to make facec planar
      par(11)=0.    
      npar=11
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)      
      chname='THSN'
      chshap='TRAP'
      nmed=med_PLASTIC_SCINT
      par(1)=1.905
      par(2)=0.0
      par(3)=0.0
      par(4)=5.795
      par(5)=6.206
      par(6)=4.211
      par(7)=0.0
      par(8)=5.795
      par(9)=7.000
      par(10)=5.005
      par(11)=0.0    
      npar=11
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)      
      chname='PMTE'
      chshap='TUBE'
      nmed=med_PYREX
      par(1)=3.7
      par(2)=3.81
      par(3)=11.43
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
c now we put all these volumes in the division CEL1
c with the right rotation matrices
c the following is the rottation matrix for LGTH relative to CEL1-first par
      LGTHR1=3000
      thetax=90.
      phix=105.
      thetay=40.
      phiy=195.
      thetaz=50.
      phiz=15.
      call GSROTM(LGTHR1,thetax,phix,thetay,phiy,thetaz,phiz)
c the same for the secnd par
      LGTHR2=LGTHR1+10
      thetax=90.
      phix=75.
      thetay=40.
      phiy=165.
      thetaz=50.
      phiz=-15.
      call GSROTM(LGTHR2,thetax,phix,thetay,phiy,thetaz,phiz)      
c rotmat for THSN and PMTE are the same      
      THSNR1=LGTHR1
      PMTER1=LGTHR1
      THSNR2=LGTHR2
      PMTER2=LGTHR2
c now positioning these elements into CEL1-first copies
      chname='THSN'
      chmoth='CEL1'      
      nr=1
      irot=THSNR1
      chonly='ONLY'
#ifdef SOLARIS
      x1=1.38*1.275*15.639*SIND(50.) 
      y1=0.0
      z=1.38*1.275*15.639*COSD(50.)-delz/2
      x=x1*COSD(15.)-y1*SIND(15.)
      y=y1*COSD(15.)+x1*SIND(15.)
#else
      x1=1.38*1.275*15.639*DSIN(50.d0) 
      y1=0.0
      z=1.38*1.275*15.639*DCOS(50.d0)-delz/2
      x=x1*DCOS(15.d0)-y1*DSIN(15.d0)
      y=y1*DCOS(15.d0)+x1*DSIN(15.d0)
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='THSN'
      chmoth='CEL1'      
      nr=2
      irot=THSNR2
      chonly='ONLY'
#ifdef SOLARIS
      x1=1.38*1.275*15.639*SIND(50.) 
      y1=0.0
      z=1.38*1.275*15.639*COSD(50.)-delz/2
      x=x1*COSD(15.)+y1*SIND(15.)
      y=y1*COSD(15.)-x1*SIND(15.)
#else
      x1=1.38*1.275*15.639*DSIN(50.d0) 
      y1=0.0
      z=1.38*1.275*15.639*DCOS(50.d0)-delz/2
      x=x1*DCOS(15.d0)+y1*DSIN(15.d0)
      y=y1*DCOS(15.d0)-x1*DSIN(15.d0)
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)      
      chname='LGTH'
      chmoth='CEL1'      
      nr=1
      irot=LGTHR1
      chonly='ONLY'
#ifdef SOLARIS
      x1=1.416*1.38*1.275*15.639*SIND(50.) ! hard wired temporarily
      y1=0.0
      z=1.416*1.38*1.275*15.639*COSD(50.)-delz/2
      x=x1*COSD(15.)-y1*SIND(15.)
      y=y1*COSD(15.)+x1*SIND(15.)
#else
      x1=1.416*1.38*1.275*15.639*DSIN(50.d0) ! hard wired temporarily
      y1=0.0
      z=1.416*1.38*1.275*15.639*DCOS(50.d0)-delz/2
      x=x1*DCOS(15.d0)-y1*DSIN(15.d0)
      y=y1*DCOS(15.d0)+x1*DSIN(15.d0)
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)      
      chname='LGTH'
      chmoth='CEL1'      
      nr=2
      irot=LGTHR2
      chonly='ONLY'
#ifdef SOLARIS
      x1=1.416*1.38*1.275*15.639*SIND(50.) ! hard wired temporarily
      y1=0.0
      z=1.416*1.38*1.275*15.639*COSD(50.)-delz/2
      x=x1*COSD(15.)+y1*SIND(15.)
      y=y1*COSD(15.)-x1*SIND(15.)            
#else
      x1=1.416*1.38*1.275*15.639*DSIN(50.d0) ! hard wired temporarily
      y1=0.0
      z=1.416*1.38*1.275*15.639*DCOS(50.d0)-delz/2
      x=x1*DCOS(15.d0)+y1*DSIN(15.d0)
      y=y1*DCOS(15.d0)-x1*DSIN(15.d0)            
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='PMTE'
      chmoth='CEL1'      
      nr=1
      irot=PMTER1
      chonly='ONLY'
#ifdef SOLARIS
      x1=2.108*1.38*1.275*15.639*SIND(50.)  ! hard wired temporarily
      y1=0.0
      z=2.108*1.38*1.275*15.639*COSD(50.)-delz/2
      x=x1*COSD(15.)-y1*SIND(15.)
      y=y1*COSD(15.)+x1*SIND(15.)                  
#else
      x1=2.108*1.38*1.275*15.639*DSIN(50.d0)  ! hard wired temporarily
      y1=0.0
      z=2.108*1.38*1.275*15.639*DCOS(50.d0)-delz/2
      x=x1*DCOS(15.d0)-y1*DSIN(15.d0)
      y=y1*DCOS(15.d0)+x1*DSIN(15.d0)                  
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      chname='PMTE'
      chmoth='CEL1'      
      nr=2
      irot=PMTER2
      chonly='ONLY'
#ifdef SOLARIS
      x1=2.108*1.38*1.275*15.639*SIND(50.)  ! hard wired temporarily
      y1=0.0
      z=2.108*1.38*1.275*15.639*COSD(50.)-delz/2
      x=x1*COSD(15.)+y1*SIND(15.)
      y=y1*COSD(15.)-x1*SIND(15.)                  
#else
      x1=2.108*1.38*1.275*15.639*DSIN(50.d0)  ! hard wired temporarily
      y1=0.0
      z=2.108*1.38*1.275*15.639*DCOS(50.d0)-delz/2
      x=x1*DCOS(15.d0)+y1*DSIN(15.d0)
      y=y1*DCOS(15.d0)-x1*DSIN(15.d0)                  
#endif
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)      
c
      end
