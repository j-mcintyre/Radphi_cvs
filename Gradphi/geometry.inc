c ** Edit History **
c
c -RTJ, Dec 15 2000
c  Craig Steffen has determined that the distance from the geometric center
c  of the target to the front surface of the lead glass wall is 102.9 cm
c  for the summer 2000 run.  To insure that I set lgdz according to
c        lgdz = 102.9 + z_shift + blockz/2 + targz + tLong/2
c  where z_shift (from placelgd.F) is 10.0cm and targz is the z position
c  of the upstream face of the target.
c

c
c these constants are material numbers
c
      integer 	mat_AIR,		mat_PLAST,
     + 		mat_CsI,		mat_ACRYLIC,
     + 		mat_CARBON,		mat_MYLAR,
     + 		mat_TUNGSTEN,		mat_ALUMINIUM,
     + 		mat_LEAD,		mat_CH2,
     + 		mat_VACUUM,		mat_IRON,
     + 		mat_PYREX,		mat_SHIELD,
     + 		mat_LEAD_GLASS_SF3,	mat_LIQUID_DEUTERIUM,
     + 		mat_BERYLLIUM,		mat_HELIUM,
     + 		mat_LIQUID_HYDROGEN,	mat_liqhelium, 		
     +		mat_liqnitrogen,	mat_liqneon,
     +		mat_copper,		mat_uranium,
     +		mat_lithium,		mat_silicon,
     +		mat_G10,		mat_nickel,
     +		mat_lead_scifi
      integer mat_carbon_fiber_epoxy, mat_scint_wrap

c           
      parameter (mat_AIR=15,		mat_PLAST=18,
     + 		mat_CsI=19,		mat_ACRYLIC=20,
     + 		mat_CARBON=6,		mat_MYLAR=22,
     + 		mat_TUNGSTEN=12,	mat_ALUMINIUM=9,
     + 		mat_LEAD=13,		mat_CH2=26,	
     + 		mat_VACUUM=16,		mat_IRON=10,
     + 		mat_PYREX=29,		mat_SHIELD=30,
     + 		mat_LEAD_GLASS_SF3=34,	mat_LIQUID_DEUTERIUM=2,
     + 		mat_BERYLLIUM=5,	mat_HELIUM=38,
     + 		mat_LIQUID_HYDROGEN=1,	mat_liqhelium=3, 		
     +		mat_liqnitrogen=7,	mat_liqneon=8,
     +		mat_copper=11,		mat_uranium=14,
     +		mat_lithium=4,		mat_silicon=46,
     +		mat_G10=48,		mat_nickel=49,
     +      mat_lead_scifi=51)
      parameter (mat_carbon_fiber_epoxy=52, mat_scint_wrap=53)
c
c these constants are media numbers
c
      integer	med_VACUUM,		med_AIR,
     + 		med_PLASTIC_SCINT,	med_ACRYLIC,
     + 		med_WAVE_L_SHI,		med_MYLAR,
     + 		med_CARBON_FIBER,	med_magnetpole,
     + 		med_TUNGSTEN,		med_IRON,
     + 		med_LEAD,		med_PYREX,
     + 		med_SHIELD,		med_BERYLLIUM_TGT,
     + 		med_ALUMINIUM,		med_LEAD_GLASS_SF3,
     + 		med_HELIUM,		med_LIQUID_HYDROGEN,
     +		med_ebeampipe,		med_gtubegas,
     +		med_gbeamtube,		med_atmosphere,
     +		med_vacwindow,		med_CLASshield,
     +		med_radiator,		med_magnetgap,
     +		med_CLAStarget, 	med_leadingap,
     +		med_leadbrick,		med_collimator,
     +      med_PbSciFi,		med_cathode,
     + 		med_wrap,			med_gap,
     +      med_CsI
      integer med_carbon_fiber_epoxy, med_scint_wrap
c     
      parameter	(med_VACUUM=1,		med_AIR=2,
     + 		med_PLASTIC_SCINT=3,	med_ACRYLIC=4,
     + 		med_WAVE_L_SHI=5,	med_MYLAR=6,
     + 		med_CARBON_FIBER=7,	med_magnetpole=37,
     + 		med_TUNGSTEN=11,	med_IRON=12,
     + 		med_LEAD=13,		med_PYREX=14,
     + 		med_SHIELD=15,		med_BERYLLIUM_TGT=22,
     + 		med_ALUMINIUM=17,	med_LEAD_GLASS_SF3=20,
     + 		med_HELIUM=23,		med_LIQUID_HYDROGEN=24,
     +		med_ebeampipe=25,	med_gtubegas=26,
     +		med_gbeamtube=27,	med_atmosphere=28,
     +		med_vacwindow=29,	med_CLASshield=30,
     +		med_radiator=31,	med_magnetgap=32,
     +		med_CLAStarget=33, 	med_leadingap=34,
     +		med_leadbrick=35,	med_collimator=36,
     +      med_PbSciFi=39,		med_cathode=40,
     + 		med_wrap=8,			med_gap=9,
     +      med_CsI=10)
      parameter (med_carbon_fiber_epoxy=41, med_scint_wrap=42)
c
c the following parameters are to position units of the detector
c namely UPV , RPD , CPV, LGD 
      real	upvx,		upvy,		upvz,
     +		upvtheta,	upvphi,		upvpsi,
     + 		rpdx,		rpdy,		rpdz,
     +		rpdtheta,	rpdphi,		rpdpsi,
     + 		cpvx,		cpvy,		cpvz,
     +		cpvtheta,	cpvphi,		cpvpsi,
     + 		lgdx,		lgdy,		lgdz,
     +		lgdtheta,	lgdphi,		lgdpsi,
     + 		targx,		targy,		targz,
     +		trgttheta,	trgtphi,	trgtpsi,
     +		wallx,		wally,		wallz     
c
      parameter	(upvx=0.,	upvy=0.,	upvz=-90.00,
     +		upvtheta=0.,	upvphi=0.,	upvpsi=0.,
     + 		rpdx=0.,	rpdy=0.,	rpdz=16.55,
     +		rpdtheta=0.,	rpdphi=15.,	rpdpsi=0.,
     + 		cpvx=0.,	cpvy=0.,	cpvz=92.0,
     +		cpvtheta=0.,	cpvphi=0.,	cpvpsi=0.,
#if defined LGD_FOR_GLUEX
     + 		lgdx=0.,	lgdy=0.,	lgdz=546.5,
#else
     + 		lgdx=0.,	lgdy=0.,	lgdz=135.3,
#endif
     +		lgdtheta=0.,	lgdphi=0.,	lgdpsi=0.,
     + 		targx=0.,	targy=0.,	targz=0.,
     +		trgttheta=0.,	trgtphi=0.,	trgtpsi=0.,
     +		wallx=0.,	wally=0.,	wallz=-111.)
c     
c these are to specify complete blocks (surrounded by kevlar)     
      real      blockx,		blocky,		blockz
      parameter	(blockx=4.0,	blocky=4.0,	blockz=45.)
c      
c these are sizes of real lg blocks
      real      blockx1,	blocky1,	blockz1
      parameter	(blockx1=3.95,	blocky1=3.95,	blockz1=45.)
c     
      integer	Nx,		Ny
#if defined LGD_FOR_GLUEX
      parameter (Nx=60,		Ny=60)
#else
      parameter (Nx=28,		Ny=28)
#endif
c
c Following are the position coordinates for herculite cover
c
      real	phfrox1,	phfroy1,	phfroz1,
     +		phsidx1,	phsidy1,	phsidz1,
     +		phtopx1,	phtopy1,	phtopz1,
     +		phfrox2,	phfroy2,	phfroz2,
     +		phsidx2,	phsidy2,	phsidz2,
     +		phtopx2,	phtopy2,	phtopz2
      parameter	(phfroy1=0.,	phfrox1=0.,	phfroz1=80.0,
     +		phsidy1=0.,	phsidx1=98.1,	phsidz1=142.4,
     +		phtopy1=116.,	phtopx1=0.,	phtopz1=142.4,
     +		phfroy2=0.,	phfrox2=0.,	phfroz2=215.1,
     +		phsidy2=0.,	phsidx2=-98.1,	phsidz2=142.4,
     +		phtopy2=-116.,	phtopx2=0.,	phtopz2=142.4)      
c
c frame around lgd detector
c dimensions
      real	frsdx,		frsdy,		frsdz,
     +		frhrx,		frhry,		frhrz,
     +		frtrx,		frtry,		frtrz
      parameter	(frsdx=10.16,	frsdy=231.,	frsdz=10.16,
     +		frhrx=10.16,	frhry=10.16,	frhrz=124.46,
     +		frtrx=175.26,	frtry=10.16,	frtrz=10.16)
c
c positions
      real	pfrsdx1,	pfrsdy1,	pfrsdz1,
     +		pfrhrx1,	pfrhry1,	pfrhrz1,
     +		pfrtrx1,	pfrtry1,	pfrtrz1,
     +		pfrsdx2,	pfrsdy2,	pfrsdz2,
     +		pfrhrx2,	pfrhry2,	pfrhrz2,
     +		pfrtrx2,	pfrtry2,	pfrtrz2,
     +		pfrsdx3,	pfrsdy3,	pfrsdz3,
     +		pfrhrx3,	pfrhry3,	pfrhrz3,
     +		pfrtrx3,	pfrtry3,	pfrtrz3,
     +		pfrsdx4,	pfrsdy4,	pfrsdz4,
     +		pfrhrx4,	pfrhry4,	pfrhrz4,
     +		pfrtrx4,	pfrtry4,	pfrtrz4
      parameter	(pfrsdy1=0.,	pfrsdx1=92.71,	pfrsdz1=75.0,
     +		pfrhry1=110.42,	pfrhrx1=92.71,	pfrhrz1=142.5,
     +		pfrtry1=110.42,	pfrtrx1=0.,	pfrtrz1=75.0,
     +		pfrsdy2=0.,	pfrsdx2=-92.71,	pfrsdz2=75.0,
     +		pfrhry2=110.42,	pfrhrx2=-92.71,	pfrhrz2=142.5,
     +		pfrtry2=-110.42,pfrtrx2=0.,	pfrtrz2=75.0,
     +		pfrsdy3=0.,	pfrsdx3=-92.71,	pfrsdz3=210.0,
     +		pfrhry3=-110.42,pfrhrx3=92.71,	pfrhrz3=142.5,
     +		pfrtry3=110.42,	pfrtrx3=0.,	pfrtrz3=210.0,
     +		pfrsdy4=0.,	pfrsdx4=92.71,	pfrsdz4=210.0,
     +		pfrhry4=-110.42,pfrhrx4=-92.71,	pfrhrz4=142.5,
     +		pfrtry4=-110.42,pfrtrx4=0.,	pfrtrz4=210.0)
c
c platform for lgd and cpv
c
      real	pla1x,		pla1y,		pla1z,
     +		pla2x,		pla2y,		pla2z 
      parameter (pla1x=175.26,	pla1y=11.43,	pla1z=0.5,
     + 		pla2x=175.26,	pla2y=0.5,	pla2z=139.830)
c      		
c parameters for placing the platform
c
      real	ppla1x,		ppla1y,		ppla1z,
     +		ppla2x,		ppla2y,		ppla2z 
#if defined LGD_FOR_GLUEX
      parameter (ppla1x=0.,	ppla1y=-132.215,ppla1z=75.25,
     + 		ppla2x=0.,	ppla2y=-126.5,	ppla2z=145.165)
#else
      parameter (ppla1x=0.,	ppla1y=-82.215,	ppla1z=75.25,
     + 		ppla2x=0.,	ppla2y=-76.5,	ppla2z=145.165)
#endif
c           
c these are parameters to define and position other volumes
c such as beampipe, sweeping magnets etc.
c      
      real	hfrox,hfroy,hfroz
      real	hsidx,hsidy,hsidz
      real	htopx,htopy,htopz
      real	eLong,eDiam,eThick
      real 	gLong,gDiam,GThick
      real	dDiam1,dDiam2,dLength
      real 	wThick
      real 	rThick,rDiam
      real 	tLong,tDiam,tDist
      real 	cLong,cDiam,cDist
      real 	hLong,hDiam
      real 	oLong,oDiam1,oDiam2,oThick
      real 	sThick,sWide,sHigh,sDiam,sDist      
      real	mWide,mHigh,mLong,mGap
      real	aLong,aBore1,aBore2
      real	aDist1,aDist2
      real	xWide,xHigh,xLong,xGap,xDiam
      real      mField,xField
      parameter (hfroy=231.14,hfrox=195.58,hfroz=0.04)	! herculite
      parameter (hsidy=231.14,hsidx=0.04,hsidz=144.78)	! cover
      parameter (htopy=0.04,htopx=195.58,htopz=144.78)	! three sheets
      parameter (eLong=22.1e2,eDiam=7.5,eThick=0.2)	! electron beam pipe
      parameter (gLong=17.6e2,gDiam=9.7,gThick=0.050)	! photon beam tube
      parameter (dDiam1=6.5,dDiam2=22.0,dLength=10.8)	! Doughnut
      parameter (wThick=100.e-4)			! vacuum windows
      parameter (rThick=0.01,rDiam=1.)			! brems radiator
      parameter (tLong=2.619,tDiam=2.873,tDist=40.5e2)	! RadPhi target
      parameter (cLong=18.5,cDiam=2.54,cDist=21.9e2)	! CLAS target
      parameter (oDiam1=7.6,oDiam2=10.2)		! CLAS shield
      parameter (oLong=150.,oThick=20.)
      parameter (sThick=10.16,sWide=100.,sHigh=100.)	! RadPhi lead wall
      parameter (sDiam=15.5,sDist=38.6e2)
      parameter (hLong=100.e2,hDiam=2.e2)		! exptl. area
      parameter (mWide=8.,mHigh=6.,mLong=20.,mGap=4.0)	! sweeping magnets
      parameter (aLong=30.,aBore1=0.8,aBore2=1.6)	! collimator appertures
      parameter (aDist1=5.80e2,aDist2=8.84e2)
      parameter (xWide=25.,xHigh=25.,xLong=50.)		! post-CLAS
      parameter (xGap=10.,xDiam=8.)			! sweeping magnet
      parameter (mField=5.,xField=8.)			! magnetic fields
          
     	
c This section of geometry.inc is for the summer 99 update of the monte carlo
c Written 6/02/99 Tom Bogue
c 
c This section contains parameters to describe the BSD
c NOTE: Some of this information is also in the file BSDgeom.h for 
c compatability w/ c routines.
c
c Inner radius for the three scintillator rings
      real ScintInRad_BSD(3)
c     data ScintInRad_BSD /35.0,36.0,37.0/
      data ScintInRad_BSD /35.3,36.0,36.7/
c outer radius for the three scintillator rings
c note the 0.5 cm gap between the rings
      real ScintOutRad_BSD(3)
c     data ScintOutRad_BSD /35.5,36.5,37.5/
      data ScintOutRad_BSD /35.8,36.5,37.2/
c half-height of the scintillators
      real ScintHeight_BSD
      data ScintHeight_BSD /35.0/
c inner/outer radius of the inner substructure wall
      real SubInWallInRad_BSD
      data SubInWallInRad_BSD /33.95/
      real SubInWallOutRad_BSD
      data SubInWallOutRad_BSD /34.05/
c inner/outer radius of the outer substructure wall
      real SubOutWallInRad_BSD
      data SubOutWallInRad_BSD /34.90/
      real SubOutWallOutRad_BSD
      data SubOutWallOutRad_BSD /35.0/
c half-height of the substructure
      real SubHeight_BSD
      data SubHeight_BSD /36.9/
c note: substructure end spans from SubInWallInRad_BSD to SubOutWallOutRad_BSD
c
c half-thickness of the end of the substructure
      real SubEndHeight_BSD
c     data SubEndHeight_BSD /3.5/
      data SubEndHeight_BSD /0.5/
c the z position of the center of the BSD in the global volume
      real Z_BSD
      data Z_BSD /31.5/
c 
c Now to describe the parameters of the BGV
c
c First the total length of a lead section
      real LeadLength_BGV
      data LeadLength_BGV /86.5/
c Now the thickness of a section
      real Thickness_BGV
      data Thickness_BGV /9.3/
c Now the length and height of the notch in each section
      real NotchHeight_BGV
      real NotchLength_BGV
      data NotchHeight_BGV /4.5/
      data NotchLength_BGV /4.5/
c Now the inner radius of the whole BGV
      real InRad_BGV
      data InRad_BGV /39.0/
c Now the thickness of the light guide in the back
      real LGthickness_BGV
      data LGthickness_BGV /1.0/
c NOTE: LGthickness is just a guess.  I'll find out for sure eventually
c and change if nesessary
c now the number of sections to the BGV
      integer NumSections_BGV
      data    NumSections_BGV /24/
c and finally the z-pos of the BGV opening
      real Z_BGV
c     data Z_BGV /-10.7/
      data Z_BGV /-7.7/		!moved 4.16.2005 [rtj]
c
c now for the parameter of the BGV frame
c
c the inner and outer radii
      real FrameInRad_BGV
      real FrameOutRad_BGV
      data FrameInRad_BGV /44.5/
      data FrameOutRad_BGV /100.0/
c the thickness of the steel and aluminium disks
      real ThicknessSteelDisk_BGV
      real ThicknessAlDisk_BGV
      data ThicknessSteelDisk_BGV /1.0/
      data ThicknessAlDisk_BGV /1.5/
c the spacing between the disks
      real DiskSpacing_BGV
      data DiskSpacing_BGV /2.0/
c The z position of the frame in global coodinates
      real FrameZ_BGV
      data FrameZ_BGV /-44.25/
c I don't know this yet

c UPV section UPV description was rewritten 6/99 by T Bogue
c
c the length, width and height of the scintillator
      real LengthScint_UPV
      parameter (LengthScint_UPV = 35.56)
      real WidthScint_UPV
      parameter (WidthScint_UPV = 17.78)
      real HeightScint_UPV
      parameter (HeightScint_UPV = 0.3175)
c the inner radius, outer radius and length of the phototube
      real InRadPMT_UPV
      parameter (InRadPMT_UPV = 3.0)
      real OutRadPMT_UPV
      parameter (OutRadPMT_UPV = 3.175)
      real LengthPMT_UPV
      parameter (LengthPMT_UPV = 17.78)
c the length of the light guide
      real LengthLG_UPV
      parameter (LengthLG_UPV = 15.24)
c the width/height of the PMT end of the lightguide
      real HeightLG_UPV
      parameter (HeightLG_UPV= 5.08)
c 
c now for the positioning of the blades w/in the UPV
c a blade is defined as a Scint, a lightguide, and a PMT
c a.k.a. VOLU 'UPVA'
c 
c the horizontal gap between the horizontal blades
      real xGapHorBlades_UPV
      parameter (xGapHorBlades_UPV = 11.43)
c the vertical gap between the vertical blades
      real yGapVertBlades_UPV
      parameter (yGapVertBlades_UPV = 8.89)
c the vertical spacing between the horizontal blades
      real ySpacingHorBlades_UPV
      parameter (ySpacingHorBlades_UPV = 15.15575)
c the z shift from center of the horizontal blades
c note that the center blades are shifted downstream,
c and the other horizontals are shifted upstream
      real zShiftHorBlades_UPV
      parameter (zShiftHorBlades_UPV = 3.81)

c CPV section, only a partial description
      real LengthScint_CPVA
      parameter (LengthScint_CPVA = 22.438*2.54)
      real LengthScint_CPVB
      parameter (LengthScint_CPVB = 17.438*2.54)
      real LengthScint_CPVC
      parameter (LengthScint_CPVC = 22.250*2.54)
      real LengthScint_CPVD
      parameter (LengthScint_CPVD = 20.625*2.54)
      real LengthScint_CPVE
      parameter (LengthScint_CPVE = 20.625*2.54)
