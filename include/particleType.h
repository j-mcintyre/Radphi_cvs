/*
 * particleType.h
*/

#ifndef particleTypeH_INCLUDED
#define particleTypeH_INCLUDED

#ifndef __CINT__
static const char sccsid_particleTypeH[] = "@(#)particleType.h\t4.1\tCreated 10 Dec 1995 18:48:51, \tcompiled "__DATE__;
#endif

typedef enum {

  /*
   * These constants are defined to be
   * same as GEANT. see http://wwwcn.cern.ch/asdoc/geant/H2GEANTCONS300.html
   * for more details.
  */

  Unknown     =  0,
  Gamma       =  1,
  Positron    =  2,
  Electron    =  3,
  Neutrino    =  4,
  MuonPlus    =  5,
  MuonMinus   =  6,
  Pi0         =  7,
  PiPlus      =  8,
  PiMinus     =  9,
  KLong       = 10,
  KPlus       = 11,
  KMinus      = 12,
  Neutron     = 13,
  Proton      = 14,
  AntiProton  = 15,
  KShort      = 16,
  Eta         = 17,
  Lambda      = 18,
  SigmaPlus   = 19,
  Sigma0      = 20,
  SigmaMinus  = 21,
  Xi0         = 22,
  XiMinus     = 23,
  OmegaMinus  = 24,
  AntiNeutron = 25,

  /* the constants defined by GEANT end here */

  /* These are E852-defined constants */

  Rho0        = 57,
  RhoPlus     = 58,
  RhoMinus    = 59,
  omega       = 60,
  EtaPrime    = 61,
  PhiMeson    = 62,

  /* These are E94-016-defined constants */

  a0_980      = 63,
  f0_980      = 64,
  
  /* sigma->Pi0Pi0 (1.0 GeV wide), Eta->3Pi0 (sigma Pi0) */
  
  sigma       = 65,
  Eta3Pi0     = 66,
 /* EtaPrime -> omega(Pi0gamma)gamma */
  EtaPrime4g  = 67
	  
} Particle_t;

#define STABLE_PARTICLE_TYPES 60

char *ParticleType(Particle_t);
float ParticleMass(Particle_t);
int   ParticleCharge(Particle_t);

#endif
/* end file */
