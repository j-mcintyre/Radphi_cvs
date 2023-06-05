c subroutine to define the UPV geometry to GEANT
c completely rewritten 6/99 by Tom Bogue
c
      subroutine upv
      implicit none
c the parameters of the volumes to be defined
      real par(12)
c a error return value for GSVOLU
c on error, the routine upvErr is called
      integer ivolu
c the position to place a volume
      real x,y,z
c the width of 1 blade
      real heightBlade
c the length of a blade
      real lengthBlade
c rotation matrixes for the left,right,top, and bottom blades
      integer forLeft,forRight,forTop,forBottom
c and here are the parameters for defining those rotation matrixes
      real thetax,phix,thetay,phiy,thetaz,phiz
c include file defining all geometry parameters for the UPV
#     include "geometry.inc"
c first off, we must define the UPV box itself
c the width of the box = x gap/2 + width of scint + lightguide + PMT
      par(1)=xGapHorBlades_UPV/2+LengthScint_UPV+LengthLG_UPV+
     +       LengthPMT_UPV
c the height of the box = y gap/2 + same
      par(2)=yGapVertBlades_UPV/2+LengthScint_UPV+LengthLG_UPV+
     +       LengthPMT_UPV
c now we need the width of 1 blade
      if (2*OutRadPMT_UPV .gt. HeightLG_UPV) then
        heightBlade=2*OutRadPMT_UPV
      else
        heightBlade=HeightLG_UPV
      endif
c the depth of the box = depth of lightguide+zshift of horizontal blades
      par(3)=heightBlade/2.0+zShiftHorBlades_UPV
c now we can define the UPV
      CALL GSVOLU('UPV ','BOX ',med_AIR,par,3,ivolu)
      if (ivolu .le. 0) call upvErr
c now to define the blade
c we already calculated the width of a blade, so
      par(1)=WidthScint_UPV/2.0
      par(2)=heightBlade/2.0
      lengthBlade=LengthScint_UPV+LengthLG_UPV+
     +             LengthPMT_UPV
      par(3)=lengthBlade/2.0
      CALL GSVOLU('UPVA','BOX ',med_AIR,par,3,ivolu)
      if (ivolu .le. 0) call upvErr
c now to define the scintillator
      par(1)=WidthScint_UPV/2.0
      par(2)=HeightScint_UPV/2.0
      par(3)=LengthScint_UPV/2.0
      CALL GSVOLU('UPVS','BOX ',med_PLASTIC_SCINT,par,3,ivolu)
      if (ivolu .le. 0) call upvErr
c now for the lightguide
      par(1)=WidthScint_UPV/2.0
      par(2)=HeightLG_UPV/2.0
      par(3)=HeightScint_UPV/2.0
      par(4)=HeightLG_UPV/2.0
      par(5)=LengthLG_UPV/2.0
      CALL GSVOLU('UPVL','TRD2',med_ACRYLIC,par,5,ivolu)
      if (ivolu .le. 0) call upvErr
c and finally for the PMT
      par(1)=InRadPMT_UPV
      par(2)=OutRadPMT_UPV
      par(3)=LengthPMT_UPV/2.0
      CALL GSVOLU('UPVP','TUBE',med_PYREX,par,3,ivolu)
      if (ivolu .le. 0) call upvErr
c now we position the scint, LG, and PMT into the blade
c first for the scint
      x=0
      y=0
      z=-(LengthPMT_UPV+LengthLG_UPV)/2.0
      CALL GSPOS('UPVS',1,'UPVA',x,y,z,0,'ONLY')
c now for the lightguide
      x=0
      y=0
      z=-(LengthPMT_UPV-LengthScint_UPV)/2.0
      CALL GSPOS('UPVL',1,'UPVA',x,y,z,0,'ONLY')
c and finally the PMT
      x=0
      y=0
      z=(LengthLG_UPV+LengthScint_UPV)/2.0
      CALL GSPOS('UPVP',1,'UPVA',x,y,z,0,'ONLY')
c now we must find the rotation matrixes to position the various blades
      forRight=8000
c the local x points along the global -y
      thetax=270
      phix=90
c the local y points along the global z
      thetay=0
      phiy=0
c the local z points along the global -x
      thetaz=90
      phiz=180
      CALL GSROTM(forRight,thetax,phix,thetay,phiy,thetaz,phiz)
c now for the next rotation
      forLeft=8100
c the local x points along the global y
      thetax=90
      phix=90
c the local y points along the global z
      thetay=0
      phiy=0
c the local z points along the global x
      thetaz=90
      phiz=0
      CALL GSROTM(forLeft,thetax,phix,thetay,phiy,thetaz,phiz)
c now for the next rotation
      forTop=8200
c the local x points along the global -x
      thetax=90
      phix=180
c the local y points along the global z
      thetay=0
      phiy=0
c the local z points along the global y
      thetaz=90
      phiz=90
      CALL GSROTM(forTop,thetax,phix,thetay,phiy,thetaz,phiz)
c now for the last rotation matrix
      forBottom=8300
c the local x points along the global x
      thetax=90
      phix=0
c the local y points along the global z
      thetay=0
      phiy=0
c the local z points along the global -y
      thetaz=90
      phiz=270
      CALL GSROTM(forBottom,thetax,phix,thetay,phiy,thetaz,phiz)
c now we can position the blades inside the UPV
c note that the copy number equals the channel number
c first the beam right blades
      x=-(xGapHorBlades_UPV+lengthBlade)/2.0
      y=-ySpacingHorBlades_UPV
      z=-zShiftHorBlades_UPV
      CALL GSPOS('UPVA',0,'UPV ',x,y,z,forRight,'ONLY')
c x doesn't change
      x=x
      y=0
      z=zShiftHorBlades_UPV
      CALL GSPOS('UPVA',2,'UPV ',x,y,z,forRight,'ONLY')
c x still doesn't change
      x=x
      y=ySpacingHorBlades_UPV
      z=-zShiftHorBlades_UPV
      CALL GSPOS('UPVA',4,'UPV ',x,y,z,forRight,'ONLY')
c now for the beam left side
      x=(xGapHorBlades_UPV+lengthBlade)/2.0
      y=-ySpacingHorBlades_UPV
      z=-zShiftHorBlades_UPV
      CALL GSPOS('UPVA',1,'UPV ',x,y,z,forLeft,'ONLY')
c again, x doesn't change
      x=x
      y=0
      z=zShiftHorBlades_UPV
      CALL GSPOS('UPVA',3,'UPV ',x,y,z,forLeft,'ONLY')
c one final time x remains the same
      x=x
      y=ySpacingHorBlades_UPV
      z=-zShiftHorBlades_UPV
      CALL GSPOS('UPVA',5,'UPV ',x,y,z,forLeft,'ONLY')
c now for the bottom blade
      x=0
      y=-(yGapVertBlades_UPV+lengthBlade)/2.0
      z=0
      CALL GSPOS('UPVA',6,'UPV ',x,y,z,forBottom,'ONLY')
c and finally the top blade
      x=0
      y=(yGapVertBlades_UPV+lengthBlade)/2.0
      z=0
      CALL GSPOS('UPVA',7,'UPV ',x,y,z,forTop,'ONLY')
      end
c
c here is a little error-handler called if GSVOLU returns ivolu .le. 0
c
      subroutine upvErr
      implicit none
      print *,'Error defining a volume in subroutine upv!'
      stop
      end

