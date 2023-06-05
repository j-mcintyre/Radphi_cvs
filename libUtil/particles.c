#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include <particleType.h>

/*
  $Log: particles.c,v $
  Revision 1.7  2004/05/28 18:18:44  radphi
  particles.c
       - roll back to the version last committed by Mihajlo [rtj]

  Revision 1.5  2004/03/24 14:16:36  radphi
  Added eta(3pi0) and etap(omega,gamma) decays - MK

  Revision 1.4  2003/04/08 16:23:49  radphi
  - added new meson decays sigma->Pi0,Pi0 and Eta3Pi0->sigma,Pi0 [mk]

  Revision 1.3  2001/03/07 00:46:27  radphi
  Changes made by jonesrt@zeus
  added several new particle types to particles.c  -rtj

 * Revision 1.2  1996/12/10  20:43:12  lfcrob
 * Added <math.h> for HUGE_VAL on non-SGI machines
 *
 * Revision 1.1  1996/11/25  20:38:01  lfcrob
 * Initial revision
 *
*/

static const char rcsid[]="$Id: particles.c,v 1.7 2004/05/28 18:18:44 radphi Exp $";



char *ParticleType(Particle_t p)
{
  static char ret[20];
  switch (p) {
  case Unknown:
    strcpy(ret,"unknown");
    break;
  case Gamma:
    strcpy(ret,"gamma");
    break;
  case Positron:
    strcpy(ret,"positron");
    break;
  case Electron:
    strcpy(ret,"electron");
    break;
  case Neutrino:
    strcpy(ret,"neutrino");
    break;
  case MuonPlus:
    strcpy(ret,"mu+");
    break;
  case MuonMinus:
    strcpy(ret,"mu-");
    break;
  case Pi0:
    strcpy(ret,"pi0");
    break;
  case PiPlus:
    strcpy(ret,"pi+");
    break;
  case PiMinus:
    strcpy(ret,"pi-");
    break;
  case KLong:
    strcpy(ret,"K0L");
    break;
  case KPlus:
    strcpy(ret,"K+");
    break;
  case KMinus:
    strcpy(ret,"K-");
    break;
  case Neutron:
    strcpy(ret,"n");
    break;
  case Proton:
    strcpy(ret,"p");
    break;
  case AntiProton:
    strcpy(ret,"pbar");
    break;
  case KShort:
    strcpy(ret,"K0S");
    break;
  case Eta:
    strcpy(ret,"eta");
    break;
  case Lambda:
    strcpy(ret,"lambda");
    break;
  case SigmaPlus:
    strcpy(ret,"Sigma+");
    break;
  case Sigma0:
    strcpy(ret,"Sigma0");
    break;
  case SigmaMinus:
    strcpy(ret,"Sigma-");
    break;
  case Xi0:
    strcpy(ret,"Xi0");
    break;
  case XiMinus:
    strcpy(ret,"Xi-");
    break;
  case OmegaMinus:
    strcpy(ret,"Omega-");
    break;
  case AntiNeutron:
    strcpy(ret,"nbar");
    break;
  case Rho0:
    strcpy(ret,"rho0");
    break;  
  case RhoPlus:
    strcpy(ret,"rho+");
    break;
  case RhoMinus:
    strcpy(ret,"rho;");
    break;
  case omega:
    strcpy(ret,"omega");
    break;
  case EtaPrime:
    strcpy(ret,"etaprime");
    break;
  case EtaPrime4g:
    strcpy(ret,"etaprime4g");
    break;
  case PhiMeson:
    strcpy(ret,"phi");
    break;
  case sigma:
    strcpy(ret,"sigma");
    break;
  case Eta3Pi0:
    strcpy(ret,"Eta3Pi0");
    break;
  case a0_980:
    strcpy(ret,"a0_980");
    break;
  case f0_980:
    strcpy(ret,"f0_980");
    break;
  default:
    sprintf(ret,"type(%d)",(int)p);
    break;
  }
  return(ret);
}

 
    
float ParticleMass(Particle_t p)
{
  switch (p) {
  default:
    fprintf(stderr,"ParticleMass: Error: Unknown particle type %d,",p);
    fprintf(stderr," returning HUGE_VAL...\n");
    return HUGE_VAL;
  case Unknown:     return HUGE_VAL;
  case Gamma:       return 0;
  case Positron:    return 0.0005101;
  case Electron:    return 0.0005101;
  case Neutrino:    return 0;
  case MuonPlus:    return 0.105658;
  case MuonMinus:   return 0.105658;
  case Pi0:         return 0.13497;
  case PiPlus:      return 0.139568;
  case PiMinus:     return 0.139568;
  case KShort:      return 0.497671;
  case KLong:       return 0.497671;
  case KPlus:       return 0.49364;
  case KMinus:      return 0.49364;
  case Neutron:     return 0.93956;
  case Proton:      return 0.93827;
  case AntiProton:  return 0.93827;
  case Eta:         return 0.54745;
  case Lambda:      return 1.11568;
  case SigmaPlus:   return 1.18937;
  case Sigma0:      return 1.19264;
  case SigmaMinus:  return 1.18937;
  case Xi0:         return 1.31483;
  case XiMinus:     return 1.32131;
  case OmegaMinus:  return 1.67245;
  case AntiNeutron: return 0.93956;
  case Rho0:        return 0.7693;
  case RhoPlus:     return 0.7693;
  case RhoMinus:    return 0.7693;
  case omega:       return 0.78257;
  case EtaPrime:    return 0.95778;
  case EtaPrime4g:  return 0.95778;
  case PhiMeson:    return 1.01942;
  case sigma:	    return 1.00000;
  case Eta3Pi0:     return 0.54745;		    
  case a0_980:      return 0.980;
  case f0_980:      return 0.980;
  }
}

int ParticleCharge(Particle_t p)
{
  switch (p) {
  default:
    fprintf(stderr,"ParticleCharge: Error: Unknown particle type %d,",p);
    fprintf(stderr," returning 0...\n");
    return 0;
  case Unknown:     return  0;
  case Gamma:       return  0;
  case Positron:    return +1;
  case Electron:    return -1;
  case Neutrino:    return  0;
  case MuonPlus:    return +1;
  case MuonMinus:   return -1;
  case Pi0:         return  0;
  case PiPlus:      return +1;
  case PiMinus:     return -1;
  case KShort:      return  0;
  case KLong:       return  0;
  case KPlus:       return +1;
  case KMinus:      return -1;
  case Neutron:     return  0;
  case Proton:      return +1;
  case AntiProton:  return -1;
  case Eta:         return  0;
  case Lambda:      return  0;
  case SigmaPlus:   return +1;
  case Sigma0:      return  0;
  case SigmaMinus:  return -1;
  case Xi0:         return  0;
  case XiMinus:     return -1;
  case OmegaMinus:  return -1;
  case AntiNeutron: return  0;
  case Rho0:        return  0;
  case RhoPlus:     return +1;
  case RhoMinus:    return -1;
  case omega:       return  0;
  case EtaPrime:    return  0;
  case EtaPrime4g:  return  0;
  case PhiMeson:    return  0;
  case sigma:       return  0;
  case Eta3Pi0:     return  0;
  case a0_980:      return  0;
  case f0_980:      return  0;
  }
}
