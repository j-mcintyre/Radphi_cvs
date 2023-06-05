c November 18, 2002
c   Introduced switch FOLLOW_CERENKOV to enable a simulation mode
c   where each Cerenkov photon is produced individually in the LGD
c   and followed using ray optics to the point of detection or
c   absorption.  To enable this feature, uncomment the following line.
c #define FOLLOW_CERENKOV
c - Richard Jones
c------------------------------------------------------------------
c Defines tracking media
c Author: Andriy Kurylov
c
      subroutine meddef
      implicit none
#     include "geometry.inc"
      integer itmed,nmat,isvol,ifield,nwbuf
      real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf
      character*20 natmed
c vacuum
      itmed=med_VACUUM
      nmat=mat_VACUUM
      natmed='vacuum'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c air
      itmed=med_AIR
      nmat=mat_AIR
      natmed='air'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c plastic scint
      itmed=med_PLASTIC_SCINT
      nmat=mat_PLAST
      natmed='plastic-scint'
      isvol=1
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c acrylic
      itmed=med_ACRYLIC
      nmat=mat_ACRYLIC
      natmed='acrylic'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c wave l. shi
      itmed=med_WAVE_L_SHI
      nmat=mat_PLAST
      natmed='wave-l.-shi'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c mylar
      itmed=med_MYLAR
      nmat=mat_MYLAR
      natmed='mylar'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c carbon fybers
      itmed=med_CARBON_FIBER
      nmat=mat_CARBON
      natmed='carbon-fiber'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c lead-glass block wrapping
      itmed=med_wrap
      nmat=mat_CH2
      natmed='wrap'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c space between lead-glass block and wrapping
      itmed=med_gap
      nmat=mat_AIR
      natmed='gap'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c CsI
      itmed=med_CsI
      nmat=mat_CsI
      natmed='CsI'
      isvol=1
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c tungsten
      itmed=med_TUNGSTEN
      nmat=mat_TUNGSTEN
      natmed='tungsten'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c iron
      itmed=med_IRON
      nmat=mat_IRON
      natmed='iron'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c lead
      itmed=med_LEAD
      nmat=mat_LEAD
      natmed='lead'
      isvol=0
      ifield=0
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c pyrex
      itmed=med_PYREX
      nmat=mat_PYREX
      natmed='pyrex'
      isvol=0
      ifield=0
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c phototube cathode
      itmed=med_CATHODE
      nmat=mat_PYREX
      natmed='photocathode'
      isvol=1
      ifield=0
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c shield - looks like concrete, mostly dirt
      itmed=med_SHIELD
      nmat=mat_SHIELD
      natmed='shield'
      isvol=0
      ifield=0
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c aluminum
      itmed=med_ALUMINIUM
      nmat=mat_ALUMINIUM
      natmed='aluminum'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)     

c lead glass SF3
      itmed=med_LEAD_GLASS_SF3
      nmat=mat_LEAD_GLASS_SF3
      natmed='lead-glass-sf3'
      isvol=1
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c beryllium target
      itmed=med_BERYLLIUM_TGT
      nmat=mat_BERYLLIUM
      natmed='beryllium'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c helium
      itmed=med_HELIUM
      nmat=mat_HELIUM
      natmed='helium'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c liquid hydrogen
      itmed=med_LIQUID_HYDROGEN
      nmat=mat_LIQUID_HYDROGEN
      natmed='liquid-hydrogen'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c electron beam pipe
      itmed=med_ebeampipe
      nmat=mat_IRON
      natmed='e-beampipe'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c photon beam pipe
      itmed=med_gbeamtube
      nmat=mat_MYLAR
      natmed='g-beamtube'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c photon beam tube content
      itmed=med_gtubegas
      nmat=mat_HELIUM
      natmed='g-tubegas'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c entrance/exit windows
      itmed=med_vacwindow
      nmat=mat_MYLAR
      natmed='vac-window'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)


c bremsstrahlung radiator
      itmed=med_radiator
      nmat=mat_CARBON
      natmed='brems-radiator'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c CLAS target
      itmed=med_CLAStarget
      nmat=mat_LIQUID_HYDROGEN
      natmed='CLAS-target'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c lead shielding material
      itmed=med_leadbrick
      nmat=mat_LEAD
      natmed='lead-brick'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c CLAS magnet shield
      itmed=med_CLASshield
      nmat=mat_IRON
      natmed='CLAS-shield'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c sweeping magnet poles
      itmed=med_magnetpole
      nmat=mat_IRON
      natmed='magnet-poles'
      isvol=0
      ifield=2
      fieldm=max(mField,xField)
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c sweeping magnet gap (field region in vacuum)
      itmed=med_magnetgap
      nmat=mat_VACUUM
      natmed='magnet-gap'
      isvol=0
      ifield=2
      fieldm=max(mField,xField)
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c lead inside sweeping magnet gap
      itmed=med_leadingap
      nmat=mat_LEAD
      natmed='lead-in-gap'
      isvol=0
      ifield=2
      fieldm=xField
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)

c photon beam collimators
      itmed=med_collimator
      nmat=mat_nickel
      natmed='collimator'
      isvol=0
      ifield=0
      fieldm=0.
      tmaxfd=-1.	! these parameters are automatically
      stemax=-1.	! set by GEANT (in its great wisdom) 
      deemax=-1.	! so don't think that setting them
      stmin=-1.		! can have any big effect.
      epsil=10.e-4
      ubuf=0.
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     &            deemax,epsil,stmin,ubuf,nwbuf)            
c
c This section by Tom Bogue summer 99
c
c medium of the lead_scifi used in the BGV
      itmed=med_PbSciFi
      nmat=mat_lead_scifi
      natmed='PbSciFi'
      isvol=1
      ifield=0
      fieldm=0
      tmaxfd=-1
      stemax=-1
      deemax=-1
      stmin=-1
      epsil=10e-4
      ubuf=0
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)

c
c begin section for summer 2020 update
c
c carbon fiber in the BSD
      itmed=med_carbon_fiber_epoxy
      nmat=mat_carbon_fiber_epoxy
      natmed='CarbonFiberEpoxy'
      isvol=0
      ifield=0
      fieldm=0
      tmaxfd=-1
      stemax=-1
      deemax=-1
      stmin=-1
      epsil=10e-4
      ubuf=0
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)
c
c wrapping around the bsd scints
      itmed=med_scint_wrap
      nmat=mat_scint_wrap
      natmed='ScintWrap'
      isvol=0
      ifield=0
      fieldm=0
      tmaxfd=-1
      stemax=-1
      deemax=-1
      stmin=-1
      epsil=10e-4
      ubuf=0
      nwbuf=0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)

#if defined FOLLOW_CERENKOV
c
c extra definitions for Cerenkov simulation -rtj
c
      call setupCerenkov
#endif
      end
