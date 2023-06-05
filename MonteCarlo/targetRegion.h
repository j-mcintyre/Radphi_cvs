/*====================================================================

Various and sundry parameters for the E-016-94 simulation.


   coordinate system:
     (0,0,0) is the center of the upstream face of the target
    
   units:
     length: cm
     energy: GeV
      angle: degrees
   momentum: GeV/c

PCF
22 JAN 96

======================================================================*/

#define GLASS_TO_TARGET_Z 100.00
#define GAMMA_E_MIN 0.150
#define THETA_MIN 40.0
#define THETA_MAX 60.0
#define MIN_GAMMA_SEP 8.00
#define MAX_CUTS 20

/*----------------------- beam */





/* ---------------------- target (cylindrical) */

const float target_radius = 1.4;
const float target_length = 2.54;

const float brem_to_target_dist = 3740.0;  /* radiator to target distance */

const float BE_RADIATION_LENGTH =  35.3; // 35.4 for Be-- 1.0 for test


/* ---------------------- pretrigger PMT */


//p=.239 corrosponds to T=.030, which should get the proton into the 
//second scintillator.  set 19JUN96 CzorS
const float minRecoil_momentum = 0.239;
const float maxRecoil_momentum = 99999.0;

const float target_to_pmt_dist = 20.0;
const float pmt_angle_degrees  = 50.0;

const int BEND = 1; 
// BEND=1, scattering; BEND=0, no scattering
const int ELOSS = 1;
// ELOSS=1, proton loses energy in target
// ELOSS=0, proton does not lose energy in target












