c A geometrical description of the BSD for the GEANT library
c written 6/2/99 by Tom Bogue
c describes the scintillator layers and the supporting structure beneath them.
c other parts should be out of the way of detected particles.
c
      subroutine bsd
      implicit none
c geometry parameters
#     include "geometry.inc"
c variables used for this routine
      real par(3)
      real Height
      integer ivolu
c the shift of the scints from center of the BSD unit
      real zShift
c PrintParams_BSD is a debuging flag that makes bsd print the parameters
c of all its volumes
c#define PrintParams_BSD
c
c first describe the outer layers of scintillator
c these are simplified to be cylindrical layers (TUBS).
c the spiraling stucture of the individual scintillators will be accounted for
c    during GUSTEP
c
      if (SubHeight_BSD .gt. ScintHeight_BSD) then
        Height=SubHeight_BSD
      else
        Height=ScintHeight_BSD
      endif
c
c the first part of this is defining the region they will occupy
      par(1)=SubInWallInRad_BSD
      par(2)=ScintOutRad_BSD(3)+0.2
      par(3)=Height
      call GSVOLU('BSD ','TUBE',med_AIR,par,3,ivolu)
c add wrapping material in between and around the scintillator layers
c august 24, 2020 -rtj-
c      call GSVOLU('BSD ','TUBE',med_scint_wrap,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
c now to define the 3 scintilator rings
      par(1)=ScintInRad_BSD(1)
      par(2)=ScintOutRad_BSD(1)
      par(3)=ScintHeight_BSD
      call GSVOLU('BSR1','TUBE',med_PLASTIC_SCINT,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
      par(1)=ScintInRad_BSD(2)
      par(2)=ScintOutRad_BSD(2)
      call GSVOLU('BSR2','TUBE',med_PLASTIC_SCINT,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
      par(1)=ScintInRad_BSD(3)
      par(2)=ScintOutRad_BSD(3)
      call GSVOLU('BSR3','TUBE',med_PLASTIC_SCINT,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
c Now to place the rings into the detector
c
      zShift=SubHeight_BSD-ScintHeight_BSD
      call GSPOS('BSR1',1,'BSD ',0,0,zShift,0,'ONLY')
      call GSPOS('BSR2',1,'BSD ',0,0,zShift,0,'ONLY')
      call GSPOS('BSR3',1,'BSD ',0,0,zShift,0,'ONLY')
c
c The next step is to describe the substructure
c which consists of a hollow cylindrical shell
c
c first define the volume they will occupy
      par(1)=SubInWallInRad_BSD
      par(2)=SubOutWallOutRad_BSD
      par(3)=SubHeight_BSD
      call GSVOLU('BSDS','TUBE',med_AIR,par,3,ivolu)
c fill the interior volume with carbon fiber instead of air
c august 24, 2020 -rtj-
c      call GSVOLU('BSDS','TUBE',med_carbon_fiber_epoxy,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
c now define the inner wall
c
      par(1)=SubInWallInRad_BSD
      par(2)=SubInWallOutRad_BSD
      par(3)=SubHeight_BSD
      call GSVOLU('BSSI','TUBE',med_ALUMINIUM,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
c Now define the outer wall
c
      par(1)=SubOutWallInRad_BSD
      par(2)=SubOutWallOutRad_BSD
      par(3)=SubHeight_BSD
      call GSVOLU('BSSO','TUBE',med_ALUMINIUM,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
c Now define the top
c
      par(1)=SubInWallInRad_BSD
      par(2)=SubOutWallOutRad_BSD
      par(3)=SubEndHeight_BSD
      call GSVOLU('BSST','TUBE',med_ALUMINIUM,par,3,ivolu)
      if (ivolu .le. 0) call BSDerr
#ifdef PrintParams_BSD
      call GPVOLU(ivolu)
#endif
c
c Note that the bottom is a copy of the top
c Now to place these parts into the substructure volume 
c
      call GSPOS('BSSI',1,'BSDS',0,0,0,0,'ONLY')
      call GSPOS('BSSO',1,'BSDS',0,0,0,0,'ONLY')
      call GSPOS('BSST',1,'BSDS',0,0,
     +           SubHeight_BSD-SubEndHeight_BSD,0,'ONLY')
      call GSPOS('BSST',2,'BSDS',0,0,
     +           -(SubHeight_BSD-SubEndHeight_BSD),0,'ONLY')
c these peices overlap but are materially identical, so it doesn't matter
c
c now to place the substructure into the BSD
      call GSPOS('BSDS',1,'BSD ',0,0,0,0,'ONLY')
      end
c
c subroutine called if there is an error in the routine BSD
c
      subroutine BSDerr
      implicit none
      print *,'Error initializing the BSD geometry in the routine bsd'
      stop
      end
