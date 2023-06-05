c This subroutine describes the geometry of the BGV to GEANT
c Written 6/03/99 by Tom Bogue
c
      subroutine bgv
      implicit none
c geometrical data is in geometry.inc
#     include "geometry.inc"
c variables used
      real par(16)
      real FullWidth
c the full width of the BGV frame
      integer ivolu
c error return variable.  GSVOLU sets it to .le. 0 if it has a problem
c if it is set to .le. 0, then the routine BGVerr is called
c and displays an error message
c PrintParams_BGV is a debuging flag that prints all volume parameters
c in bgv if defined
c#define PrintParams_BGV
c
c The first step is to define the volume of the detector
c
      par(1)=-180.0/real(NumSections_BGV)
      par(2)=360.0
      par(3)=NumSections_BGV
      par(4)=2
c
      par(5)=0
      par(6)=InRad_BGV
      par(7)=InRad_BGV+Thickness_BGV
c
      par(8)=LeadLength_BGV+LGthickness_BGV
      par(9)=par(6)
      par(10)=par(7)
      call GSVOLU('BGV ','PGON',med_PbSciFi,par,10,ivolu)
      if (ivolu .le. 0) call BGVerr
#ifdef PrintParams_BGV
      call GPVOLU(ivolu)
#endif
c
c Now to divide the volume into individual sections
      call GSDVN('BGVS','BGV ',NumSections_BGV,2)
c
c Now to define the light guide on the back of the lead sections
c
      par(2)=360.0/NumSections_BGV
      par(1)=-par(2)/2.0
      par(3)=1
      par(4)=4
c
      par(5)=LeadLength_BGV-NotchLength_BGV
      par(6)=InRad_BGV
      par(7)=par(6)
c
      par(8)=LeadLength_BGV
      par(9)=par(6)
      par(10)=par(7)+NotchHeight_BGV
c
      par(11)=par(8)
      par(12)=par(9)
      par(13)=InRad_BGV+Thickness_BGV
c
      par(14)=par(11)+LGthickness_BGV
      par(15)=par(12)
      par(16)=par(13)
      call GSVOLU('BVLG','PGON',med_ACRYLIC,par,16,ivolu)
      if (ivolu .le. 0) call BGVerr
#ifdef PrintParams_BGV
      call GPVOLU(ivolu)
#endif
c Now to place the light guide inside of the lead section
      call GSPOS('BVLG',1,'BGVS',0,0,0,0,'ONLY')
c
c Now on to the frame
c
c First find the full width of the frame
      FullWidth=ThicknessSteelDisk_BGV+ThicknessAlDisk_BGV+
     +          DiskSpacing_BGV
c now define the frame volume
c
      par(1)=FrameInRad_BGV
      par(2)=FrameOutRad_BGV
      par(3)=FullWidth/2.0
      call GSVOLU('BGVF','TUBE',med_air,par,3,ivolu)
      if (ivolu .le. 0) call BGVerr
#ifdef PrintParams_BGV
      call GPVOLU(ivolu)
#endif
c Now define the steel disk
      par(3)=ThicknessSteelDisk_BGV/2.0
      call GSVOLU('BVFS','TUBE',med_IRON,par,3,ivolu)
      if (ivolu .le. 0) call BGVerr
#ifdef PrintParams_BGV
      call GPVOLU(ivolu)
#endif
c now define the Al disk
      par(3)=ThicknessAlDisk_BGV/2.0
      call GSVOLU('BVFA','TUBE',med_ALUMINIUM,par,3,ivolu)
      if (ivolu .le. 0) call BGVerr
#ifdef PrintParams_BGV
      call GPVOLU(ivolu)
#endif
c now place the Al disk into the frame volume
      call GSPOS('BVFA',1,'BGVF',0,0,
     +           (FullWidth-ThicknessAlDisk_BGV)/2.0,0,'ONLY')
c now for the Steel disk
      call GSPOS('BVFS',1,'BGVF',0,0,
     +           (FullWidth-ThicknessSteelDisk_BGV)/2.0,0,'ONLY')
c now to put the frame into the global system
c     call GSPOS('BGVF',1,'HALL',0,0,FrameZ_BGV,0,'ONLY')
      end
c
      subroutine BGVerr
      implicit none
c called if GSVOLU sets ivolu .le. 0
      print *,'Error in BGV volume declaration is subroutine BGV'
      stop
      end
