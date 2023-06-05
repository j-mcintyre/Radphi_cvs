/* The event display here represented was designed, coded & documented by   */
/* Troels Petersen under the supervision of Elton Smith, with help from     */
/*   Dave  Heddle and James Muguira.					    */
/* Comments can be found around more complicated code leaving the basics    */
/* uncommented. Upon problems and/or questions, please contact Troels       */
/* Petersen at: <petersen@cebaf.gov> or <petersen@fys.ku.dk>		    */

#include <Hv.h>
#include <itypes.h>

/* - initialize the global variables for the drawing package - */

int	pedestals=1;
int	mam_display=1;

Hv_Item	LeadGlassItem, TaggerItem, RecoilItem,
	BGVItem, CPVItem, BSDItem, UPVItem; 
Hv_Item	Box1, Box2, Box3, Box4, Box5, Box6;

/* - the size (in pixels) of each LGB - */
int	nsize=12;
/* - pixels per centimeter ( = (float)nsize/4.0 ) - */
float	ppc=3.0;
/* - the detector geometry; geo[2]=7 means that 2nd row has 2*7 blocks. (see 'init.c') - */
int	geo[15];
/* - the current position of the mouse (in cm), (see feedback.c) - */
float	worldx, worldy;

/* - energy, momenta and position of the photons hitting the detector - */
float	e[8], px[8], py[8], pz[8], x[8], y[8], z[8];

/* - same as above, for Monte Carlo generated values - */
float	eMC[8], pxMC[8], pyMC[8], pzMC[8], xMC[8], yMC[8], zMC[8];

/* - number of each kind of meson reconstructed - */
int 	npi, neta, nrho, nomega, netap, nphi;

/*  the number of marked photons and the current event number - */
int	nmarked=0;
int16   trigbit=0;
int16   trigflag=0;
int16   trigval=0;
int16   onlineflag=0;
Boolean calibrated=True;

Boolean matchup=False;
Boolean lgscale=False;
Boolean phodata=False;

typedef struct {
  int multicount;
  float time[16];
} Tlist;

int     run=0;
int     event=0;
int     mam_energy=0;
int     mam_mass=0;
float	lg[28][28];
float	rpd_g[6], rpd_e[12];
float 	bsd_r[12], bsd_l[12], bsd_s[24];
Tlist	bsd_r_times[12], bsd_l_times[12], bsd_s_times[28];
float 	bgv_ds[24], bgv_us[24];
Tlist 	bgv_ds_times[24], bgv_us_times[24];
float 	cpv[30];
Tlist 	cpv_times[30];
float 	upv[8];
Tlist 	upv_times[9];
Tlist 	tagger_times[39];

int nmarked, mark[8], nhit, noMC;
float etphot, miphot, tiphot, ebeam;

/* the energy of beam in GeV default value: 5.65 GeV */
float   beamEnergy = 5.65;

/* - struct for multihit TDC info - */

typedef struct {
  int times[16];
  int multicount;
} tdc_stuff;

/* - mark[] tells whether or not a photon has been marked by the mouse - */
int	mark[8];

/* - sets bounds on valid T and E for each detector - */
/* - T in ns, E in GeV, deltaE in mips - */
float lb_cpvT = -10.,	ub_cpvT = 15.;
float lb_upvT = -10.,	ub_upvT = 15.;
float lb_bgvT = -10.,	ub_bgvT = 30.;
float lb_bsdT = -10.,	ub_bsdT = 30.;
float lb_tagT = -10.,	ub_tagT = 10.;
float lb_upvE = 1.,	ub_upvE = 10.;
float lb_cpvE = 1.,	ub_cpvE = 10.;
float lb_bsdE = 1.,	ub_bsdE = 10.;
float lb_bgvE = 1.,	ub_bgvE = 10.;
float lb_lgdE = 1.,	ub_lgdE = 10.;

/* - x- & y-position of the center and the upper left corner of the detector (in pixels) - */
short	xc, yc, xp, yp;
short   xctcpv, yctcpv, xctag, yctag, xct, yct, xctupv, yctupv, xctbsd, yctbsd;
/* - the name of the event path+file, the filename alone and the filepointer to it - */
char	*filename, *fname = "null";

float epi,mpi,tpi,p1pi,p2pi,eta,p1eta,p2eta,Omega,p2omega,phi,p2phi;
