#ifndef   __EDH
#define   __EDH

/* adjustable event display parameters */

#define LGD_RED_HOT	400
#define BSD_RED_HOT	400
#define BGV_RED_HOT	400
#define CPV_RED_HOT	400
#define UPV_RED_HOT	400
#define TAGGER_FULL_SCALE	400

/* Set bounds on valid TDC,ADC ranges */

#define LB_ADC 1
#define UB_ADC 4096
#define LB_BSDTDC 1
#define UB_BSDTDC 2200
#define LB_BGVTDC 1
#define UB_BGVTDC 2200
#define LB_CPVTDC 1
#define UB_CPVTDC 2200
#define LB_UPVTDC 1
#define UB_UPVTDC 2200

/* define view tags */

#define  ED_VIEW     1
#define  TAGGER_VIEW 2
#define  CPV_VIEW 3
#define  UPV_VIEW 4
#define  BSD_VIEW 5

/* define feedback tags */

#define WORLDX        0
#define WORLDY        1
#define DIST          2
#define THETA         3
#define POSS          4
#define INFO          5

typedef struct {
  int multicount;
  float time[16];
} Tlist;

#include <Hv.h>
#include <ntypes.h>
#include <mctypes.h>
#include <disData.h>

extern Hv_Widget optped, optnoped;  
extern Hv_Item	Box1, Box2, Box3, Box4, Box5, Box6;
extern Hv_Item	LeadGlassItem, TaggerItem, RecoilItem,
                BGVItem, CPVItem, BSDItem, UPVItem; 

extern float	lb_cpvT, ub_cpvT;
extern float 	lb_upvT, ub_upvT;
extern float	lb_bgvT, ub_bgvT;
extern float	lb_bsdT, ub_bsdT;
extern float	lb_tagT, ub_tagT;
extern float	lb_cpvE, ub_cpvE;
extern float	lb_upvE, ub_upvE;
extern float	lb_bgvE, ub_bgvE;
extern float	lb_bsdE, ub_bsdE;
extern float	lb_lgdE, ub_lgdE;

/* the number of pixels used per block in LGD display */
extern int	nsize;
/* - pixels per centimeter ( = (float)nsize/4.0 ) - */
extern float    ppc;
/* the detector geometry; geo[2]=7 means that 2nd row has 2*7 blocks. 
 *   (see 'init.c')                                                   */
extern int      geo[15];
/* the current position of the mouse (in cm), (see feedback.c) */
extern float    worldx, worldy;
/* coordinates used for pointer feedback */
extern short    xc, yc, xp, yp, xctag, yctag;
extern short    xctcpv, yctcpv, xctbsd, yctbsd, xctupv, yctupv;


/* matchup, lgscale & phodata determines whether the routines are
 * used in the drawing
 */
extern Boolean  matchup, lgscale, phodata;

/* switch for calibrated energies and trigger selector */
extern short 	trigbit;
extern short 	trigflag;
extern short 	trigval;
extern short	onlineflag;
extern Boolean	calibrated;
extern int	run, event;

/* - energy, momenta and position of the photons hitting the detector - */
extern float	e[8], px[8], py[8], pz[8], x[8], y[8], z[8];

/* - same as above, for Monte Carlo generated values - */
extern float	eMC[8], pxMC[8], pyMC[8], pzMC[8], xMC[8], yMC[8], zMC[8];

/* - the number of each type of meson reconstructed - */
/* - while 0 < npi < 4, the number of other specific mesons should not exceed 1 - */
extern int 	npi, neta, nrho, nomega, netap, nphi;

/* - energy sum, inv. mass and inv. t of reconstructed Pi and other mesons - */
extern float	epi[4], mpi[4], tpi[4], eta[3], rho[3], Omega[3],
		etap[3], phi[3]; 

/* - photon number of reconstructed mesons - */
extern int	p1pi[4], p2pi[4], p1eta, p2eta, p1rho, p2rho, p1omega,
		p2omega, p1etap, p2etap, p1phi, p2phi;

/* fileopen tells whether a file has been opened or not */
extern Boolean  fileopen;

/* flag to signal monte carlo data (based on presence of GROUP_MC_EVENT) */
extern int monte_carlo;

extern float beamEnergy;           /* the energy of beam in GeV */
extern float lg[28][28];           /* E values for the LGD blocks */
extern float bsd_r[12],            /* deltaE values for BSD counters right, */
             bsd_l[12],            /*                                left,  */
             bsd_s[24];            /*                          and straight */
extern Tlist bsd_r_times[12],      /* T values for the BSD counters right,  */
             bsd_l_times[12],      /*                               left,   */
             bsd_s_times[28];      /*                          and straight */
extern float bgv_ds[24],           /* E values for the BGV downstream end   */
             bgv_us[24];           /*                          and upstream */
extern Tlist bgv_ds_times[24],     /* T values for the BGV downstream end   */
             bgv_us_times[24];     /*                          and upstream */
extern float cpv[30];              /* deltaE values for the CPV counters */
extern Tlist cpv_times[30];        /* T values for the CPV counters */
extern float upv[8];               /* deltaE values for the UPV counters */
extern Tlist upv_times[9];         /* T values for the UPV counters */
extern Tlist tagger_times[39];     /* T values for the tagger */
extern int   rpd_g[6], rpd_e[12];  /* adc values for RPD counters */
extern int   mam_energy, mam_mass; /* values for MAM registers */
extern int   mam_display;
extern int   pedestals;

extern int nmarked;  /* number of photons that have been marked */
extern int mark[8];  /* flags if a photon has been marked by the mouse - */
extern int nhit;     /* the number of photons hitting the detector */
extern int noMC;     /* same as above, for Monte Carlo generator */
extern float etphot, /* total energy of marked photons */
             miphot, /* invariant mass  "         "    */
             tiphot, /* t-invariant     "         "    */
             ebeam;  /* incident beam energy */

extern char *filename, /* the name of the event path+file */
            *fname;    /* the filename alone */

/* functions in init.c */
void  Init(void);

/* functions in setup.c */
void  SetupEdView(Hv_View view);
void  SetupTagView(Hv_View view);
void  SetupCPVView(Hv_View view);
void  SetupBSDView(Hv_View view);
void  SetupUPVView(Hv_View view);

/* functions in feedback.c */
void  Feedback(Hv_View, Hv_Point);
void  FeedbackTagger(Hv_View, Hv_Point);
void  FeedbackCPV(Hv_View, Hv_Point);
void  FeedbackBSD(Hv_View, Hv_Point);
void  FeedbackUPV(Hv_View, Hv_Point);

/* functions in input.c */
void  MAMSelect(Hv_Event);
void  ReadForward(Hv_Event);
void  NewFile(Hv_Event);
void  ReadBackward(Hv_Event);
void  ReadGoto(Hv_Event);
void  NewInput(Hv_Event);
void  TrigSelect(Hv_Event);
void  BeamEnergySelect(Hv_Event);
void  Online(Hv_Event);
void  EditColor(Hv_Event);
void  GetData(itape_header_t*, uint32*);
int   mcEventPrint(mc_event_t*);
int   mcVertexPrint(int, mc_vertex_t*);

/* functions in draw.c */
void  MakeNew(void);
void  LeadGlass(Hv_Item, Hv_Region);
void  Recoil(Hv_Item, Hv_Region);
void  BGV(Hv_Item, Hv_Region);
void  BSD(Hv_Item, Hv_Region);
void  Tagger(Hv_Item, Hv_Region);
void  CPV(Hv_Item, Hv_Region);
void  UPV(Hv_Item, Hv_Region);
void  BSDUnfolded(Hv_Item, Hv_Region);
void  DrawMatchup(Hv_Item, Hv_Region);
void  DrawLGscale(Hv_Item, Hv_Region);
void  DrawConnect(Hv_Item, Hv_Region);
void  DrawCPV(float,float,float,float,short);
void  DrawUPV(float,float,float,float,short);
void  PatternFill(int,int,int,short,short); 
void  Photon(Hv_Item, Hv_Region);
void  Marked(Hv_Item, Hv_Region);
void  MarkPhoton(Hv_Event);

#endif
