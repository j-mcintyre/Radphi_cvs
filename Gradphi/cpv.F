c ---------------------------------------------------------------------------
c              Charge Particle Veto
c----------------------------------------------------------------------------
c
      subroutine cpv
      implicit none
#     include "geometry.inc"
      character*4 chname,chshap
      integer nmed,npar,ivolu
      real par(12)
      character*4 chmoth,chonly
      integer nr,irot
c angles for rotational matrices      
      real thetax, phix, thetay, phiy, thetaz, phiz
c other paremeters
      real x,y,z
c these are rotation matrices      
      integer skewBR, skewBL, flipBR, flipBL
c the half-width of the box
      real Half_Width
      parameter (Half_Width=10.4)
c
c CPV 
      chname='CPV ' 
      chshap='BOX '
      nmed=med_AIR
      par(1)=172./2
      par(2)=120./2
      par(3)=Half_Width
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)      
c                  
c           
c scintillators
c
      chname='CPSA' 
      chshap='BOX '
      nmed=med_PLASTIC_SCINT
      par(1)=3.500*1.27
      par(2)=0.2
      par(3)=LengthScint_CPVA/2.0
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPSB' 
      chshap='BOX '
      nmed=med_PLASTIC_SCINT
      par(1)=6.750*1.27
      par(2)=0.2
      par(3)=LengthScint_CPVB/2.0
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPSC' 
      chshap='BOX '
      nmed=med_PLASTIC_SCINT
      par(1)=1.250*1.27
      par(2)=0.2
      par(3)=LengthScint_CPVC/2.0
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPSD' 
      chshap='BOX '
      nmed=med_PLASTIC_SCINT
      par(1)=1.250*1.27
      par(2)=0.2
      par(3)=LengthScint_CPVD/2.0
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPSE' 
      chshap='BOX '
      nmed=med_PLASTIC_SCINT
      par(1)=1.250*1.27
      par(2)=0.2
      par(3)=LengthScint_CPVE/2.0
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)

c now we define the lightguides      
      chname='CPGA' 
      chshap='TRD2'
      nmed=med_ACRYLIC
      par(1)=3.500*1.27
      par(2)=2.000*1.27
      par(3)=0.250*1.27
      par(4)=0.250*1.27
      par(5)=2.980*1.27
      npar=5
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPGB' 
      chshap='TRD2'
      nmed=med_ACRYLIC
      par(1)=6.750*1.27
      par(2)=2.000*1.27
      par(3)=0.250*1.27
      par(4)=0.250*1.27
      par(5)=6.980*1.27
      npar=5
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPGC' 
      chshap='TUBE'
      nmed=med_ACRYLIC
      par(1)=0.0
      par(2)=1.875*1.27
      par(3)=2.980*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPGD' 
      chshap='TUBE'
      nmed=med_ACRYLIC
      par(1)=0.0
      par(2)=1.875*1.27
      par(3)=3.043*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chname='CPGE' 
      chshap='TUBE'
      nmed=med_ACRYLIC
      par(1)=0.0
      par(2)=1.875*1.27
      par(3)=3.043*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)

C now CPV PMTs
      chname='CPMT' 
      chshap='TUBE'
      nmed=med_PYREX
      par(1)=1.850*1.27
      par(2)=2.000*1.27
      par(3)=5.500*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)

c the following is the PMT base      
      chname='CPMB' 
      chshap='TUBE'
      nmed=med_ACRYLIC
      par(1)=0.0
      par(2)=2.125*1.27
      par(3)=2.375*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)

c the paddles need to go diagonally into their mother to avoid overlaps
      irot=6010
      thetax=90.
      phix=315.
      thetay=90.
      phiy=45.
      thetaz=0.
      phiz=0.
      call GSROTM(irot,thetax,phix,thetay,phiy,thetaz,phiz)            

c now stick scintillators, guides, tubes and bases together in a mother
      chname='CPVA'
      chshap='BOX '
      nmed=med_AIR
      par(1)=(3.500+0.250)*1.27/sqrt(2.)
      par(2)=par(1)
      par(3)=33.042*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='CPVA'
      chname='CPSA'
      chonly='ONLY'
      nr=1
      x=0
      y=0
      z=-par(3)+(22.438)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=1
      chname='CPGA'
      z=z+(22.438+2.980)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=1
      chname='CPMT'
      z=z+(2.980+5.500)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=1
      chname='CPMB'
      z=z+(5.500+2.375)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      chname='CPVB'
      chshap='BOX '
      nmed=med_AIR
      par(1)=(6.750+0.250)*1.27/sqrt(2.)
      par(2)=par(1)
      par(3)=32.043*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='CPVB'
      chname='CPSB'
      chonly='ONLY'
      nr=1
      x=0
      y=0
      z=-par(3)+(17.438)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=1
      chname='CPGB'
      z=z+(17.438+6.980)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,irot,chonly)
      nr=2
      chname='CPMT'
      z=z+(6.980+5.500)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=2
      chname='CPMB'
      z=z+(5.500+2.375)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      chname='CPVC'
      chshap='BOX '
      nmed=med_AIR
      par(1)=2.200*1.27
      par(2)=2.200*1.27
      par(3)=33.043*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='CPVC'
      chname='CPSC'
      chonly='ONLY'
      nr=1
      x=0
      y=0
      z=-par(3)+(22.250)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=1
      chname='CPGC'
      z=z+(22.250+2.980)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=3
      chname='CPMT'
      z=z+(2.980+5.500)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=3
      chname='CPMB'
      z=z+(5.500+2.375)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      chname='CPVD'
      chshap='BOX '
      nmed=med_AIR
      par(1)=2.200*1.27
      par(2)=2.200*1.27
      par(3)=31.480*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='CPVD'
      chname='CPSD'
      chonly='ONLY'
      nr=1
      x=0
      y=0
      z=-par(3)+(20.625)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=1
      chname='CPGD'
      z=z+(20.625+3.043)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=4
      chname='CPMT'
      z=z+(3.043+5.500)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=4
      chname='CPMB'
      z=z+(5.500+2.375)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      chname='CPVE'
      chshap='BOX '
      nmed=med_AIR
      par(1)=2.200*1.27
      par(2)=2.200*1.27
      par(3)=31.480*1.27
      npar=3
      call GSVOLU(chname,chshap,nmed,par,npar,ivolu)
      chmoth='CPVE'
      chname='CPSE'
      chonly='ONLY'
      nr=1
      x=0
      y=0
      z=-par(3)+(20.625)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=1
      chname='CPGE'
      z=z+(20.625+3.043)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=5
      chname='CPMT'
      z=z+(3.043+5.500)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)
      nr=5
      chname='CPMB'
      z=z+(5.500+2.375)*1.27
      call GSPOS(chname,nr,chmoth,x,y,z,0,chonly)

c now we can put our volumes in CPV 
c because the shape is rather complicated have to position separately
c define the rotation matrices
      skewBR=6000
      thetax=45.
      phix=270.
      thetay=45.
      phiy=90.
      thetaz=90.
      phiz=180.
      call GSROTM(skewBR,thetax,phix,thetay,phiy,thetaz,phiz)
      skewBL=7000
      thetax=135.
      phix=270.
      thetay=135.
      phiy=90.
      thetaz=90.
      phiz=0.
      call GSROTM(skewBL,thetax,phix,thetay,phiy,thetaz,phiz)            
      flipBR=6100
      thetax=90.
      phix=90.
      thetay=180.
      phiy=0.
      thetaz=90.
      phiz=180.
      call GSROTM(flipBR,thetax,phix,thetay,phiy,thetaz,phiz)
      flipBL=7100
      thetax=90.
      phix=90.
      thetay=0.
      phiy=0.
      thetaz=90.
      phiz=0.
      call GSROTM(flipBL,thetax,phix,thetay,phiy,thetaz,phiz)            

c place the counters from outside in
      chmoth='CPV '
      chname='CPVB'
      chonly='ONLY'
      x=32.043*1.27
      y=19.020*2.54
      z=3.825*2.54-Half_Width
      nr=29
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      y=-y
      nr=1
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      x=-x
      nr=0
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      y=-y
      nr=28
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      chname='CPVA'
      chonly='ONLY'
      x=33.043*1.27
      y=14.026*2.54
      z=3.460*2.54-Half_Width
      nr=27
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      y=-y
      nr=3
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      x=-x
      nr=2
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      y=-y
      nr=26
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      x=33.043*1.27
      y=10.626*2.54
      z=3.095*2.54-Half_Width
      nr=25
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      y=-y
      nr=5
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      x=-x
      nr=4
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      y=-y
      nr=24
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      x=33.043*1.27
      y=7.227*2.54
      z=2.730*2.54-Half_Width
      nr=23
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      y=-y
      nr=7
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      x=-x
      nr=6
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      y=-y
      nr=22
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      x=33.043*1.27
      y=3.827*2.54
      z=2.365*2.54-Half_Width
      nr=21
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      y=-y
      nr=9
      call GSPOS(chname,nr,chmoth,x,y,z,skewBL,chonly)
      x=-x
      nr=8
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      y=-y
      nr=20
      call GSPOS(chname,nr,chmoth,x,y,z,skewBR,chonly)
      chname='CPVC'
      x=33.043*1.27
      y=2.126*2.54
      z=6.802*2.54-Half_Width
      nr=19
      call GSPOS(chname,nr,chmoth,x,y,z,flipBL,chonly)
      y=-y
      nr=11
      call GSPOS(chname,nr,chmoth,x,y,z,flipBL,chonly)
      x=-x
      nr=10
      call GSPOS(chname,nr,chmoth,x,y,z,flipBR,chonly)
      y=-y
      nr=18
      call GSPOS(chname,nr,chmoth,x,y,z,flipBR,chonly)
      chname='CPVD'
      x=(2*33.043-31.480)*1.27
      y=1.150*2.54
      z=4.428*2.54-Half_Width
      nr=17
      call GSPOS(chname,nr,chmoth,x,y,z,flipBL,chonly)
      y=-y
      nr=13
      call GSPOS(chname,nr,chmoth,x,y,z,flipBL,chonly)
      x=-x
      nr=12
      call GSPOS(chname,nr,chmoth,x,y,z,flipBR,chonly)
      y=-y
      nr=16
      call GSPOS(chname,nr,chmoth,x,y,z,flipBR,chonly)
      chname='CPVE'
      x=(2*33.043-31.480)*1.27
      y=0
      z=2.053*2.54-Half_Width
      nr=15
      call GSPOS(chname,nr,chmoth,x,y,z,flipBL,chonly)
      x=-x
      nr=14
      call GSPOS(chname,nr,chmoth,x,y,z,flipBR,chonly)
      end
