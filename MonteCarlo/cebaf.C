////////////////////////////////////////////////////////////////////
//
// cebaf.C
//
//
// Target region simulation for E-94-016.
//
//
// The geometry and acceptance code should be placed
// in CHEES once we have decided on a setup.
//
//
// COORDINATE SYSTEM:
//  0,0,0 is the center of the upstream face of the target
//
// PCF
// 30 NOV 1995
//
////////////////////////////////////////////////////////////////////

/* #define HELICITY_AD */

#include <Chees.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ntypes.h>
#include <Four.h>
#include <Three.h>
#include <SqMatrix3.h>
#include <string.h>
extern"C"{
#include <lgdGeom.h>
#include <wabbit.h>
#include <map_manager.h>
#include <dataIO.h>
#include <disData.h>
#include <tracking.h>
#include <ieeefp.h>
}

#include <util.h>
#include <targetRegion.h>
#define BUFF_SIZE 300000
#define NUM_TRACKS 20
#define MAX_GAMMA 8

#define ENERGY_SMEAR_FLOOR .02
#define ENERGY_SMEAR_STAT  .08
#define LGD_X_SMEAR_VALUE .5
#define LGD_Y_SMEAR_VALUE .5



itape_header_t *ih=NULL;
FILE *buffer,*tabulator=NULL;
extern double costheta_between(Four vec1, Four vec2);

/* some global variables */
int run_number=0;
char target_type[10];
char beam_type[10];
char recoil_type[10];
char product_type[10];
int type_set=0;
int config=0;
int total_accepted,partial_cuts_accepted,reject_reason_totals[MAX_CUTS];
int g_nAccepted = 0;                 /* counter for accepted events */
int g_nTotal = 0;                    /* counter for total events that go through target */
int g_nMissingTarget = 0;
int g_nRecoilEnLow = 0;
int g_nRecoilEnHi = 0;
int gEscapedTarget = 0;

float t_reac_true;

Four g_in_glass[NUM_TRACKS];

FILE *glassfile=NULL;


cheesInfo_t 	gFourVectors;           /* four vecs are stored here per event */
mc_param_t gParameters;        	/* general input gParameters */
particles_db_t *gParticleDBPtr;     	/* particle specifications */

float gElectronBeamE;  /* electron beam energy */

char *reasons[MAX_CUTS];

/* some function protypes */

float dot(Four vec1,Four vec2);
void doStuff(void);
void printInfo(void);
void printHelp(void);
void InitThemHistos(void);
void SaveThemHistos(void);
int  DoCebafAcceptance(double);

Three FindExitPoint(Three,Three);
int PointInTarget(Three *);
Three FindVertex(Three, Three);

int glass_sim(Three vertx,int recoil_passed);
Four smear_4v(Four invec);
Four smearPosition_4v(Four invec);
////////////////////////////////////////////////////////////////////////


void InitThemHistos(void)
{
  // make two histos for each item: accepted and generated events
  for (short i = 0; i <= 1; i++)
    {
      short prefix = 10000; // accepted plots are id 10000 + histo #
      
      /*      bbook1(i*prefix + 10,"px, beam (GeV/c) GEN",100,-0.03,0.03,0.0); // beam momenta
      bbook1(i*prefix + 11,"py, beam (GeV/c) GEN",100,-0.03,0.03,0.0);
      bbook1(i*prefix + 12,"pz, beam (GeV/c) GEN",100,0.0,4.1,0.0);
      bbook1(i*prefix + 13,"beam theta (m_e/e_max) GEN",100,0.0,6.0,0.0);
      // angle with z-axis in units of (mass electron)/electron_beam_energy

      */
      bbook2(i*prefix + 101,"p recoil (GeV/c) vs angle(deg) GEN",100,0.0,1.0,100,0.0,90.0,0.0);
      
      bbook1(i*prefix + 102,"pf, target (GeV/c) GEN",100,0.0,0.4,0.0); // mag target p
      bbook1(i*prefix + 103,"px, target (GeV/c) GEN",100,-0.3,0.3,0.0);// components
      bbook1(i*prefix + 104,"py, target (GeV/c) GEN",100,-0.3,0.3,0.0);
      bbook1(i*prefix + 105,"pz, target (GeV/c) GEN",100,-0.3,0.3,0.0);
      
      bbook1(i*prefix + 106,"recoil angle (deg) GEN",100,0.0,90.0,0.0); // angle w/ z axis
      
      /*      bbook1(i*prefix + 107,"vertex x GEN",100,0.0,target_radius,0.0);
      bbook1(i*prefix + 108,"vertex y GEN",100,0.0,target_radius,0.0);
      bbook1(i*prefix + 109,"vertex z GEN",100,0.0,target_length,0.0);
      */
      bbook1(i*prefix + 110,"recoil energy loss (GeV) GEN",100,0.0,0.3,0.0);
      bbook1(i*prefix + 111,"recoil angle (deg) mult scat GEN",100,0.0,90.0,0.0); // angle w/ z axis
      bbook2(i*prefix + 112,"p recoil (GeV/c) vs angle(deg) MS-EL GEN",100,0.0,1.0,100,0.0,90.0,0.0);
      
      bbook1(i*prefix + 113,"mult scatter theta (rad)",100,0.0,0.2,0.0); 
      
      /*      bbook1(i*prefix + 150,"recoil KE GEN",100,0.0,0.6,0.0);
      bbook1(i*prefix + 151,"recoil KE GEN, MS-EL",100,0.0,0.6,0.0);
      
      bbook2(i*prefix + 152,"KE recoil (GeV) vs angle(deg) GEN",100,0.0,1.0,100,0.0,90.0,0.0);
      bbook2(i*prefix + 153,"KE recoil (GeV) vs angle(deg) MS-EL GEN",100,0.0,1.0,100,0.0,90.0,0.0);
      
      bbook1(i*prefix + 154,"dist of recoil along PMT plane (cm)",100,0.0,40.0,0.0);
      
      bbook2(i*prefix + 155,"recoil position in PMT plane (cm, r vs z)",200,0.0,31.0,200,0.0,31.0,0.0);
      
      bbook1(i*prefix+52,"angle of impact on recoil detector cone (from beam direction) (deg)",180,0.0,90.0,0.0);
      bbook1(i*prefix+53,"p_z of incoming photon",250,0.0,4.5,0.0);
      bbook1(i*prefix+54,"mass of vec3+vec4",200,0.0,2.0,0.0);
      bbook1(i*prefix+55,"mass of vec3+vec5",200,0.0,2.0,0.0);
      bbook1(i*prefix+56,"mass of vec4+vec5",200,0.0,2.0,0.0);
      bbook1(i*prefix+57,"mass of vec3+vec4+vec5",200,0.0,2.0,0.0);
      
      
      bbook1(i*prefix+60,"mag of vector 0",200,0.0,2.0,0.0);
      bbook1(i*prefix+61,"mag of vector 1",200,0.0,2.0,0.0);
      bbook1(i*prefix+62,"mag of vector 2",200,0.0,2.0,0.0);
      bbook1(i*prefix+63,"mag of vector 3",200,0.0,2.0,0.0);
      bbook1(i*prefix+64,"mag of vector 4",200,0.0,2.0,0.0);
      bbook1(i*prefix+65,"mag of vector 5",200,0.0,2.0,0.0);
      bbook1(i*prefix+66,"mag of vector 6",200,0.0,2.0,0.0);
      */
    } /* end loop */
  
  bbook2(10058,"gamma hits on glass",100,-56.0,56.0,100,-56.0,56.0,0.0);
  bbook1(10059,"rho of gamma hits on glass",200,0.0,85.0,0.0);
  bbook2(10068,"accepted gamma hits on glass",100,-56.0,56.0,100,-56.0,56.0,0.0);
  bbook1(10069,"rho of accepted gamma hits on glass",200,0.0,85.0,0.0);
  bbook1(10067,"Energy of glass impacting gamma (GeV)",200,0.0,4.0,0.0);

  // for testing
  bbook1(50,"recoil angle from chees(deg)",200,0.0,100.0,0.0);
  bbook1(10070,"Min gamma sep in glass(cm)",150,0.0,50.0,0.0);
  bbook1(10080,"Min gamma sep in glass(cm) for accepted",150,0.0,50.0,0.0);
  bbook1(10071,"Min gamma r (cm), intermediate cuts",200,0.0,40.0,0.0);
  bbook1(10072,"Max gamma r (cm), intermediate cuts",200,0.0,65.0,0.0);
  bbook1(10081,"Min gamma r (full accepted)(cm)",200,0.0,40.0,0.0);
  bbook1(10082,"Max gamma r (full accepted)(cm)",200,0.0,65.0,0.0);

  /*  bbook1(10090,"Beam to product angle(deg)",180,0.0,180.0,0.0);
  bbook1(10091,"Product to daughter angle(deg)",180,0.0,180.0,0.0);
  bbook1(10092,"Beam to scattered photon angle(deg)",180,0.0,180.0,0.0);
  bbook1(10093,"target p to recoil p  angle(deg)",180,0.0,180.0,0.0);
  */

  /*  bbook1(10073,"Min gamma r(cm),raw",200,0.0,40.0,0.0); 
  bbook1(10074,"Max gamma r(cm),raw",200,0.0,65.0,0.0);
  bbook1(10083,"Angle of daughter with product; product rest frame(deg)",180,0.0,180.0,0.0); */
  bbook1(10084,"generated -t (GeV/c)^2",200,0.0,2.0,0.0);
  bbook1(10085,"accepted -t (GeV/c)^2",200,0.0,2.0,0.0);
  /*  bbook2(10086,"accepted PMT cone hits",100,0.0,31.0,100,0.0,31.0,0.0); */
  bbook1(10086,"generated ptsquared (GeV/c)^2",200,0.0,2.0,0.0);
  bbook1(10087,"accepted ptsquared (GeV/c)^2",200,0.0,2.0,0.0);
  
  bbook1(10088,"flight time of recoil proton (ns)",150,0.0,3.0,0.0);
  bbook1(10089,"flight time of recoil proton (accepted) (ns)",150,0.0,3.0,0.0);
  
  bbook1(10094,"accepted(NR) -t (GeV/c)^2",200,0.0,2.0,0.0);

  /*  bbook1(11054,"mass of vec3+vec4(accepted)",200,0.0,2.0,0.0);
  bbook1(11055,"mass of vec3+vec5(accepted)",200,0.0,2.0,0.0);
  bbook1(11056,"mass of vec4+vec5(accepted)",200,0.0,2.0,0.0);
  bbook1(11057,"mass of vec3+vec4+vec5(accepted)",200,0.0,2.0,0.0); */

  bbook1(12001,"-t from recoil proton (GeV/c)^2",200,0.0,2.0,0.0);
  bbook1(12003,"-t from photons (GeV/c)^2",200,0.0,2.0,0.0);

  bbook1(12004,"Energy of photons that missed glass",200,0.0,2.0,0.0);
  bbook2(12005,"Energy of photons vs. theta",100,0.0,3.2,150,0.0,4.0,0.0);

  bbook1(12010,"dN/d(costheta) (costheta) of daughter from product (gen)",100,-1.0,1.0,0.0);
  bbook1(12011,"dN/d(costheta) (costheta) of daughter from product (acc)",100,-1.0,1.0,0.0);
  
  bbook1(12012,"dN/d(costheta) (costheta) of daughter from product (recon)",100,-1.0,1.0,0.0);

  bbook2(12100,"generated pt squared vs. -t",200,0.0,1.0,200,0.0,1.0,0.0);
  bbook2(12101,"accepted pt squared vs. -t",200,0.0,1.0,200,0.0,1.0,0.0);
  
  bbook1(12200,"3 gamma generated product mass, GeV",200,0.0,2.0,0.0);
  bbook1(12201,"3 gamma smeared generated product mass, GeV",200,0.0,2.0,0.0);
  bbook1(12202,"3 gamma smeared accepted product mass, GeV",200,0.0,2.0,0.0);
  
  bbook2(12300,"accepted reconstructed product E vs. true beam E (GeV)",
	 200,0.0,6.0,200,0.0,6.0,0.0);
  bbook2(12310,"accepted reconstructed product Hel costheta vs. true",
	 100,-1.0,1.0,100,-1.0,1.0,0.0);

}


////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

float angle_4v(Four vec1,Four vec2)
{
  double one_dot_two,angle,arg;
  one_dot_two=(vec1.vec.space.x*vec2.vec.space.x)+(vec1.vec.space.y*vec2.vec.space.y)+(vec1.vec.space.z*vec2.vec.space.z);
  arg = one_dot_two/(vec1.mag3()*vec2.mag3());
  if(fabs(arg)>1){
    if (arg<0)
      return(PI);
    else
      return(0);
  }
  angle=acos(one_dot_two/(vec1.mag3()*vec2.mag3()));
  return angle;
}

Four into_rest_frame_of(Four vec1,Four vec2)
{
  Four beta,transformed;
  double energy,betavalsq,bX,bY,bZ,gamma,gm_one;
  //  double temp=vec2.vec.t;
  beta=vec2*(1/vec2.vec.t);
  bX=beta.vec.space.x;
  bY=beta.vec.space.y;
  bZ=beta.vec.space.z;
  betavalsq=SQR(beta.mag3());
  gamma=vec2.vec.t*(1/vec2.mag4());
  gm_one=gamma-1.0E+0;
  transformed.vec.t=
    (gamma*vec1.vec.t)+
    (-(gamma*bX)*vec1.vec.space.x)+
    (-(gamma*bY)*vec1.vec.space.y)+
    (-(gamma*bZ)*vec1.vec.space.z);
  transformed.vec.space.x=
    (-(gamma*bX)*vec1.vec.t)+
    ((1+(gm_one*(SQR(bX)/betavalsq)))*vec1.vec.space.x)+
    (((gm_one*bX*bY)/betavalsq)*vec1.vec.space.y)+
    (((gm_one*bX*bZ)/betavalsq)*vec1.vec.space.z);
  transformed.vec.space.y=
    (-(gamma*bY)*vec1.vec.t)+
    (((gm_one*bX*bY)/betavalsq)*vec1.vec.space.x)+
    ((1+(gm_one*(SQR(bY)/betavalsq)))*vec1.vec.space.y)+
    (((gm_one*bY*bZ)/betavalsq)*vec1.vec.space.z);
  transformed.vec.space.z=
    (-(gamma*bZ)*vec1.vec.t)+
    (((gm_one*bX*bZ)/betavalsq)*vec1.vec.space.x)+
    (((gm_one*bY*bZ)/betavalsq)*vec1.vec.space.y)+
    ((1+(gm_one*(SQR(bZ)/betavalsq)))*vec1.vec.space.z);  
  return transformed;
  
}


vector4_t convd4to4(dvector4_t invec_d)
{
  vector4_t converted;
  converted.space.x=(float)invec_d.space.x;
  converted.space.y=(float)invec_d.space.y;
  converted.space.z=(float)invec_d.space.z;
  converted.t=(float)invec_d.t;
  return converted;
}




vector3_t conv4to3(vector4_t invec)
{
  vector3_t outvec;
  outvec.x=invec.space.x;
  outvec.y=invec.space.y;
  outvec.z=invec.space.z;
  return outvec;
}  


double invmass(Four A)
{
  double tmp;
  tmp = (A.vec.t*A.vec.t)-(A.vec.space.x*A.vec.space.x)-(A.vec.space.y*A.vec.space.y)-(A.vec.space.z*A.vec.space.z);

// Hack - asssumes that any value < 0 is due to roundoff error.

  if(tmp<0.0)
    return(0);
  return sqrt(tmp);
}

void doStuff(void)   /// This is the main program
{
  float qtheta_rad,qtheta_deg,minus_t_recoil,minus_t_gammas;
  int jones,reject_flag=500;
  double x,y,z,theta,rho,dq,dqq;
  double pt_squared,smeared_product_mass,product_mass;
  double smeared_product_E,true_beam_E,reconstructed_phi_in_Hel;;
  Four recoil_in_phi,recoil,constructed_product;
  Four vec4,daughter_Ft,daughter_in_phi,beam_typical;
  Four smeared_g_in_glass[NUM_TRACKS],smeared_product;
  Four dummyFour,reconstructed_CM,b_r_CM,b_r_bat;
  static geo_tracks_rec_t *reaction_tracks=NULL;
  
  if(config==0)
    {
      printf("configuration flag not given, exiting\n\n");
      exit(5);
    }

  if(!reaction_tracks)
    reaction_tracks=(geo_tracks_rec_t*)malloc(sizeof(geo_tracks_rec_t)+((NUM_TRACKS-1)*sizeof(geo_track_t)));

  (*reaction_tracks).ntracks=NUM_TRACKS;

  //   for testing, CzorS
  
  
  ////
  //// Have CHEES generate the four-vecs for this reaction.
  //// It fills in the "gFourVectors" struct.
  ////
  //// At this point in time (JAN 96) CHEES does not do
  //// CEBAF acceptance internally-- we'll do it from this
  //// program later on.  
  
  int passed_angular_selection=0;
  int rejected;
  Four product;
  while(passed_angular_selection==0){
    rejected = GenerateOne(&gFourVectors);
    
    
    // int data_newItape(void* itape);
    data_newItape(ih);
    ih->runNo=-(run_number);
    
    if(type_set==0)
      {
	type_set=1;
	strcpy(beam_type,ParticleType(gFourVectors.data[0].id));
	strcpy(target_type,ParticleType(gFourVectors.data[1].id));
	strcpy(product_type,ParticleType(gFourVectors.data[2].id));
	/*	if(config==2)
	  strcpy(recoil_type,ParticleType(gFourVectors.data[5].id));
	else
	  if(config==3)
	    strcpy(recoil_type,ParticleType(gFourVectors.data[6].id));
	  else
	    {
	      printf("Config screwed up in typeset");
	      exit(26);
	    }
	    */
	strcpy(recoil_type,ParticleType(gFourVectors.data[3+config].id));
      }
    
    //fill in writing routines to output  12JUL96
    
    rejected = 0; // for now we say the event is accepted.  more on this later.
    
    ////
    //// Now put these four-vecs in convenient form (Four objects)
    //// so we can do things like +  and - with them
    //// 
    
    /*	 Four recoil(gFourVectors.data[5].p);  // puts recoil into a Four (four vector) object */
    
    //config
    /*    if(config==2)
      recoil=gFourVectors.data[5].p;  */
      recoil=gFourVectors.data[3+config].p;
      /*    else 
      if(config==3)
	recoil=gFourVectors.data[6].p;
      else
	printf("config screwed up for recoil assignment"); */
    product=gFourVectors.data[2].p;
    
    /* Now do Craig's kluged daughter angular distribution selector */
    Four batch_in_prod;
    Four recoil_in_prod;
    double cos_theta;
    double test_val;
    if(config==3){
#ifdef HELICITY_AD
      batch_in_prod=gFourVectors.data[5].p; /* product's batchellor photon */
      recoil_in_prod=into_rest_frame_of(recoil,product);
      batch_in_prod=into_rest_frame_of(batch_in_prod,product);
      cos_theta=costheta_between(batch_in_prod,recoil_in_prod);
      test_val=unfRand();
      test_val=test_val*2.0; /* scale to function max */
      if(test_val<(1.0+(cos_theta*cos_theta))) 
#endif
	passed_angular_selection=1;
    }
    else{
      passed_angular_selection=1;
    }
  }  /* while(passed_angular_selection==0) */
  /* ********************************************************
     ********************************************************
     
     At this point, the event is considered "official"

     ********************************************************
     ********************************************************
   */
  // puts recoil into a Four (four vector) object
  //this subscript should be 6 for phi, 3 for rho 
  //there's just _got_ to be a better way.  CzorS 20JUN96
  
  Four beam(gFourVectors.data[0].p);    // puts beam into a Four
  Four target(gFourVectors.data[1].p);  // get the target, so we know pf
  
  t_reac_true=(recoil+(target*(-1.0))).mag4Sq();
  bf1(10084,(-(t_reac_true)),1.0);
  pt_squared=SQR(product.vec.space.x)+SQR(product.vec.space.y);
  bf1(10086,pt_squared,1.0);
  bf2(12100,(-(t_reac_true)),pt_squared,1.0);

  true_beam_E=beam.vec.t;

  //Check this really well...CzorS 11JUL96
  int qz;
  /*  g_in_glass[0]=Four(gFourVectors.data[3].p);
  g_in_glass[1]=Four(gFourVectors.data[4].p);
  if(config==3){
    g_in_glass[2]=Four(gFourVectors.data[5].p);
    constructed_product=g_in_glass[0]+
                        g_in_glass[1]+
                        g_in_glass[2];
  } */

  g_in_glass[0]=Four(gFourVectors.data[3].p);
  for(qz=1;qz<config;qz++){
    g_in_glass[qz]=Four(gFourVectors.data[3+qz].p);
    /*    g_in_glass[qz].print(); */
  }

  constructed_product=g_in_glass[0];
  for(qz=1;qz<config;qz++){
    constructed_product=constructed_product+g_in_glass[qz];
  }


  bf1(12200,constructed_product.mag4(),1.0);

  //config

  smeared_product=Four(0.0,0.0,0.0,0.0);


              /*  if(config==2)
                  for(jones=0;jones<6;jones++)
                    {
                  bob[jones]=Four(gFourVectors.data[jones].p);
                  //      bf1(10060+jones,invmass(bob[jones]),1.0);
                }
              else
                    if(config==3){
                  for(jones=0;jones<7;jones++){
        	{
        	  bob[jones]=Four(gFourVectors.data[jones].p);
        	}
              } */
        	      /*  for(qz=0;qz<config;qz++){
        	     smeared_g_in_glass[qz]=smear_4v(g_in_glass[qz]);
        	     smeared_product=smeared_product+smeared_g_in_glass[qz];
        	   } */
      
                   /*    }   */
  smeared_g_in_glass[0]=smear_4v(g_in_glass[0]);
  smeared_product=smeared_g_in_glass[0];
  for(qz=1;qz<config;qz++){
    smeared_g_in_glass[qz]=smear_4v(g_in_glass[qz]);
    smeared_product=smeared_product+smeared_g_in_glass[qz];
  }
  /*  if(config==3) */
         	    /*    daughter_Ft=bob[3]+bob[4]; */

  daughter_Ft=g_in_glass[0];
  for(qz=1;qz<(config-1);qz++){
    daughter_Ft=daughter_Ft+g_in_glass[qz];
  }

  Four stationary_proton=Four(0.0,0.0,0.0,.938);
  
  smeared_product_mass=smeared_product.mag4();
  smeared_product_E=smeared_product.vec.t;

  Four reconstructed_beam=Four(0.0,0.0,smeared_product_E,smeared_product_E);

  reconstructed_CM=reconstructed_beam+stationary_proton;

  

  if(!isnand(smeared_product_mass))
    bf1(12201,smeared_product_mass,1.0);
  
  //  bf1(10054,invmass(daughter_Ft),1.0);

                             /*  if(config==3)
                              reaction_tracks->ntracks=8;
                              else 
                                if(config==2)
                                  reaction_tracks->ntracks=7; */
  /* I hope this works.  Craig Steffen 20FEB99 */
  reaction_tracks->ntracks=config+5;


  /* The following is, I believe, all to prepare for output of an 
     event to itape 
     for, for analysis by another program.  Craig Steffen 20FEB99 */

  reaction_tracks->track[0].p=conv4to3(gFourVectors.data[0].p);
  reaction_tracks->track[0].particleType=gFourVectors.data[0].id;
  reaction_tracks->track[0].charge=ParticleCharge(reaction_tracks->track[0].particleType);  
  reaction_tracks->track[0].vertex=0;
  
  reaction_tracks->track[1].p=conv4to3(gFourVectors.data[1].p);
  reaction_tracks->track[1].particleType=gFourVectors.data[1].id;
  reaction_tracks->track[1].charge=ParticleCharge(reaction_tracks->track[1].particleType);  
  reaction_tracks->track[1].vertex=0;
  
  reaction_tracks->track[2].p=conv4to3(gFourVectors.data[2].p);
  reaction_tracks->track[2].particleType=gFourVectors.data[2].id;
  //  reaction_tracks->track[2].charge=ParticleCharge(reaction_tracks->track[2].particleType);  
  reaction_tracks->track[2].vertex=0;
  
  //  reaction_tracks->track[3].p=conv4to3(convd4to4(daughter_Ft.vec));
  //  reaction_tracks->track[3].particleType=gFourVectors.data[0].id;
  //  reaction_tracks->track[3].charge=ParticleCharge(reaction_tracks->track[0].particleType);  
  //  reaction_tracks->track[3].vertex=0;
  
  reaction_tracks->track[4].p=conv4to3(gFourVectors.data[3].p);
  reaction_tracks->track[4].particleType=gFourVectors.data[3].id;
  reaction_tracks->track[4].charge=ParticleCharge(reaction_tracks->track[4].particleType);  
  reaction_tracks->track[4].vertex=0;
  
  reaction_tracks->track[5].p=conv4to3(gFourVectors.data[4].p);
  reaction_tracks->track[5].particleType=gFourVectors.data[4].id;
  reaction_tracks->track[5].charge=ParticleCharge(reaction_tracks->track[5].particleType);  
  reaction_tracks->track[5].vertex=0;
  
  //config
  reaction_tracks->track[6].p=conv4to3(gFourVectors.data[5].p);
  reaction_tracks->track[6].particleType=gFourVectors.data[5].id;
  reaction_tracks->track[6].charge=ParticleCharge(reaction_tracks->track[6].particleType);  
  reaction_tracks->track[6].vertex=0;
  
  //config CzorS
  reaction_tracks->track[7].p=conv4to3(gFourVectors.data[6].p);
  reaction_tracks->track[7].particleType=gFourVectors.data[6].id;
  reaction_tracks->track[7].charge=ParticleCharge(reaction_tracks->track[7].particleType);  
  reaction_tracks->track[7].vertex=0;
  
  /* Daughter angular distribution check */

  daughter_in_phi=into_rest_frame_of(daughter_Ft,product);
  recoil_in_phi=into_rest_frame_of(recoil,product);

  qtheta_rad=angle_4v(recoil_in_phi,daughter_in_phi);
  qtheta_deg=qtheta_rad*(180.0/PI);

  //  bf1(10083,qtheta_deg,1.0);
  bf1(12010,cos(qtheta_rad),1.0);
  //  bf1(10090,qtheta_deg,1.0);

  /*  qtheta_rad=angle_4v(bob[0],bob[5]);
  qtheta_deg=qtheta_rad*(180.0/PI); */
  //  bf1(10092,qtheta_deg,1.0);

  //  qtheta=angle_4v(bob[1],bob[6]);
  //  qtheta=qtheta*(180.0/PI);
  //  bf1(10093,qtheta,1.0);

  /*    qtheta_rad=angle_4v(product,(bob[3]+bob[4]));
    qtheta_deg=qtheta_rad*(180.0/PI); */
  //    bf1(10091,qtheta_deg,1.0);
 

  
  // testing CzorS
  vec4.vec=recoil.vec;
  x=recoil.vec.space.x;
  y=recoil.vec.space.y;
  z=recoil.vec.space.z;
  
  
  if(z==0.0)
    theta=90.0;
  else if(z>0.0) //if theta is in quadrant I
    {
      rho=sqrt(SQR(x)+SQR(y));
      theta=(atan(rho/z)*57.295);
    }
  else if(z<0.0) //if theta is in quadrant III
    {
      rho=sqrt(SQR(x)+SQR(y));
      theta=(atan(rho/z)*57.295)+180.0;
      //this is to put angle in quadrant II rather than IV
    }
  bf1(50,theta,1.0);
  
  // testing CzorS 03JUN96
  
  
  ////	    
  //// Find recoil energy/angle.
  ////	  
  
  float p_recoil = recoil.mag3();
  Four zAxis4(0.0,0.0,1.0,0.0);
  float cosTheta = (recoil.dot3(&zAxis4)) / (recoil.mag3()); 
  
  double angle = acos(cosTheta) * 360.0 / 6.28318; 
  
  ////
  //// Determine beam polar angle, in units of characteristic brem angle
  ////
  
  Three beamP(beam.vec.space);
  Three zAxis(0.0,0.0,1.0);
  
  double beamCosTheta = (zAxis * beamP)/(beamP.Mag()); // note cool c++ dot product
  double beamTheta = acos(beamCosTheta);
  double thetaC = 0.000511 / gElectronBeamE;     // mass of e- / electron beam max energy
  // is characteristic brem angle.
  beamTheta /= thetaC;
  
  ////
  //// Find the reaction vertex
  ////
  
  // first find where the beam hits the upstream end of the target
  
  double s = brem_to_target_dist / beam.vec.space.z;
  
  Three hitPoint( beam.vec.space.x * s, beam.vec.space.y * s, 0.0);
  
  // That's x,y,z coords //
  
  
  // Now choose a vertex.  First make sure hitPoint is within the target.               
  
  double r = sqrt( SQR(hitPoint.vec.x) + SQR(hitPoint.vec.y) );
  
  if (r > target_radius)
    {
      rejected = 1;
      g_nMissingTarget++;  // exit if beam missed target
      return;
    }
  
  Three beam3momentum(beam.vec.space); // make new 3-vec w/ beam momentum
  
  Three vertex = FindVertex(hitPoint,beam3momentum);
  
  ////
  //// Now model multiple scattering.
  //// we have a vertex and a recoil momentum;
  ////  determine what the exit point would be,
  ////  and then find the halfway point to this exit point.
  ////  From there we will change the direction of flight,
  ////  and find the real exit point.
  ////
  //// In pictures, we are modelling the scattering like this:
  ////
  ////    \  original path
  ////     \ 
  ////      \
      ////       \
      ////        = = = = = = =  new path	   (path bent in only one place.)
  ////
  
  Three recoilMomentum(recoil.vec.space);	
  
  Three wouldBeExit = FindExitPoint(vertex,recoilMomentum);
  
  s = (wouldBeExit.vec.z - vertex.vec.z)/ recoilMomentum.vec.z;
  
  Three bendPoint = vertex + recoilMomentum * (s/2.0);  // halfway point
  
  // Now determine new momentum
  // Rotate old momentum by theta, phi;
  //   phi uniform;
  //   theta (space angle) according to 1992 RPP (Phys Rev D) p1253.
  //               "The Particle People...Serving You Every Day!" 
  
  
  // for this model, the distance travelled through the target
  // is computed from the original (straight) path.
  
  double length = ((wouldBeExit - vertex).Mag()) / BE_RADIATION_LENGTH;
  if(length==0.0)
    length=1.0e-6;
  
  
  
  
  double beta = recoilMomentum.Mag() / recoil.vec.t; // beta = p/E
  
  double sigma_theta = (0.0136 / (beta * recoilMomentum.Mag())) * sqrt(length) *
    (1.0 + 0.038 * log(length) ) * (1.41421356237);
  
  double ms_theta = GetGaussianDist(0.0,sigma_theta);             
  double ms_phi = 2*PI*unfRand();
  
  // compute sines and cosines used in the construction of rotation
  // matrices.  See pg 41 of my CEBAF notes.
  //
  // PCF 26 JAN 96
  
  double bx,by,bz;
  bx = recoilMomentum.vec.x;
  by = recoilMomentum.vec.y;  // short-hand for formulas below.
  bz = recoilMomentum.vec.z;
  
  double costhetaXZsq = SQR(bz) / (SQR(bx) + SQR(bz)); // actually cos^2        
  double sinthetaXZsq = 1.0 - costhetaXZsq;
  double   stxz = sqrt(sinthetaXZsq);
  double   ctxz = sqrt(costhetaXZsq);
  
  // this was moved here 04JUN96 	     
  stxz = ( bx > 0.0 ? 
	   ABS(stxz) :
	   -ABS(stxz) );
  //end of moved CzorS
  
  
  
  double zprime = (bx * stxz) + (bz * ctxz);
  
  // zprime is new beam z momentum after rotation
  // about the y axis.
  
  double costhetaYZsq = SQR(zprime) / (SQR(by) + SQR(zprime)); // actually cos^2T
  double sinthetaYZsq = 1.0 - costhetaYZsq;
  double ctyz = sqrt(costhetaYZsq);
  double styz = sqrt(sinthetaYZsq);
  
  
  // Get the signs right.
  // Positive angle rotations go from the x or y axis towards the
  // z axis.  The sign makes a difference in the sine terms.
  
  
  // 	    stxz = ( bx > 0.0 ? 
  //       ABS(stxz) :
  //      -ABS(stxz) );
  
  // 	    styz = ( by > 0.0 ? 
  // ABS(styz) :
  // -ABS(styz) );
  
  //Modeified by CzorS 04JUN96
  //		  stxz = ( bx > 0.0 ? 
  //   ABS(stxz) :
  //   -ABS(stxz) );
  
  
  
  styz = ( by > 0.0 ? 
	   -ABS(styz) :
	   ABS(styz) );
  
  
  // construct matrix to align momentum w/ z axis
  // see pg. 27 of my CEBAF notes
  // PCF
  
  
  //Three CMGJmatR1(ctxz,			0.0,			-stxz);
  //Three CMGJmatR2(-styz*stxz,	ctyz,			-styz*ctxz);
  //Three CMGJmatR3(stxz*ctyz,     styz,			ctyz*ctxz);
  //CMGJmat has been modified by CzorS 04JUN96
  
  Three CMGJmatR1(ctxz,	0.0,	 -stxz);
  Three CMGJmatR2(styz*stxz,	ctyz,	 styz*ctxz);
  Three CMGJmatR3(stxz*ctyz,  -styz,	 ctyz*ctxz);
  
  
  
  SqMatrix3 CMGJMat(CMGJmatR1,CMGJmatR2,CMGJmatR3);
  
  // construct inverse of the matrix which aligns momentum
  // with z-axis
  
  //	 Three GJCMmatR1(ctxz,			-stxz*styz,		stxz*ctyz);
  //Three GJCMmatR2(0.0,			ctyz,			styz);
  //Three GJCMmatR3(-stxz,			-styz*ctxz,		ctxz*ctyz);
  
  // change made 04JUN96 CzorS
  
  Three GJCMmatR1(ctxz,			stxz*styz,		stxz*ctyz);
  Three GJCMmatR2(0.0,			ctyz,			-styz);
  Three GJCMmatR3(-stxz,			styz*ctxz,		ctxz*ctyz);
  
  
  SqMatrix3 GJCMMat(GJCMmatR1,GJCMmatR2,GJCMmatR3);
  
  
  // construct matrix to do the standard (theta,phi)
  // polar coordinate rotation
  
  double cth = cos(ms_theta);
  double sth = sin(ms_theta);
  double cfi = cos(ms_phi);
  double sfi = sin(ms_phi);
  
  Three row1( cfi*cth, -sfi, cfi*sth  );
  Three row2( sfi*cth,  cfi,  sfi*sth );
  Three row3( -sth,     0.0,      cth );
  
  SqMatrix3 BendingMatrix(row1,row2,row3);
  
  // finally, we can rotate the vector.
  // Steps: (1) align z-axis with momentum
  //        (2) rotate vector by (theta,phi) polars
  //        (3) move axes back
  
  // 	 Three newRecoilMomentum = 
  // GJCMMat * (BendingMatrix * (CMGJMat * recoilMomentum)); 
  // modified 18JUN96 CzorS
  
  Three newRecoilMomentum;
  if(BEND==1)
    {
      newRecoilMomentum = 
	GJCMMat * (BendingMatrix * (CMGJMat * recoilMomentum)); 
    }
  else if(BEND==0)
    {newRecoilMomentum = recoilMomentum;}
  else 
    {
      printf("BEND incorrectly initialized, %i\n",BEND);
      exit(-3);
    };
  
  
  //test section
  //	 Three testmom;
  // testmom=CMGJMat * recoilMomentum;
  //testmom.Print();
  //end test section CzorS
  
  
  
  // now model energy loss.
  // dE/dx = - 1.848 * 1.515 * beta ^ (-5/3)   MeV/cm
	 // see 1994 RPP p1251
  
  // remember that beta can change ==> dE/dX not constant
  
  Three newExit = FindExitPoint(bendPoint,newRecoilMomentum);
  
  double newDistance = (newExit - bendPoint).Mag() + (bendPoint - vertex).Mag();	 
  double recoilMass = recoil.mag4();	
  double currentRecoilEnergy = recoil.vec.t;
  double currentRecoilMomentum = newRecoilMomentum.Mag();
  
  double integralStepSize = newDistance / 20.0;
  // do the integral
  
  double energyChange;
  energyChange = 0;
  
  double dEdx = 0;
  if(ELOSS==1)
    {
      for ( int i = 1; i <= 20; i++)
	{  
	  beta = currentRecoilMomentum / currentRecoilEnergy;// beta = p/E
	  
	  if (beta == 0.0) break; // ran out of steam
	  
	  dEdx = (-0.001) * 1.848 * 1.515 * pow(beta,-1.6666666666);
	  
	  energyChange = dEdx * integralStepSize;
	  
	  currentRecoilEnergy += energyChange;       
	  
	  currentRecoilMomentum = SQR(currentRecoilEnergy) - SQR(recoilMass);
	  
	  if (currentRecoilMomentum < 0.0) 
	    { 
	      currentRecoilMomentum = 0.0;
	      currentRecoilEnergy = recoilMass;
	    };
	  
	  currentRecoilMomentum = sqrt(currentRecoilMomentum);
	}
    }
  else if(ELOSS==0)
    energyChange=0;
  else 
    {
      printf("ELOSS incorreclty initialized,%d\n",ELOSS);
      exit(-2);
    }

  if(beta>0.001)  //arbitrary; this is really slow, though
    bf1(10088,(.66667/beta),1.0);
  

  ///////////// 
  // -t distributions calculated here CzorS 30SEP96
  
  minus_t_recoil=2.0*recoilMass*(currentRecoilEnergy-recoilMass);
  
  //  double E_typ=.84606*gElectronBeamE;
  //  beam_typical=Four(E_typ,0.0e0,0.0e0,E_typ);
  /*  beam_typical=Four(0.0e0,0.0e0,bob[0].vec.t,bob[0].vec.t); */
  beam_typical=Four(0.0e0,0.0e0,beam.vec.t,beam.vec.t);
  
  Four subtractor;
  if(config==2)
    {
      subtractor=(g_in_glass[0]+g_in_glass[1])*(-1.0);
      minus_t_gammas=((beam_typical+subtractor).mag4Sq())*(-1.0);
    }
  else if(config==3)
    {
      subtractor=(g_in_glass[0]+g_in_glass[1]+g_in_glass[2])*(-1.0);
      minus_t_gammas=((beam_typical+subtractor).mag4Sq())*(-1.0);
    }
  
  // Now "currentRecoilMomentum" holds the recoil momentum after energy loss.       
  
  if(currentRecoilMomentum > 0.0)
    {
      gEscapedTarget = 1; // remember if we made it out or not
      // 1 = made it, 0 = died
    }
  else
    {
      gEscapedTarget = 0;
    };
  
  
  ////
  //// Now plot recoil distribution along PMT plane
  //// (so we can figure out where to put the scintillators)
  //// 
  //// see pgs 45-46 of my CEBAF notes  --------PCF
  ////
  //// Here's what's happening:
  ////  The PMT plane  is a line in r-z space.
  ////  First find the equation of this line.
  ////  Then, we use the recoil's momentum to find
  ////  the point of intersection with the PMT plane.
  ////  We plot the distribution of recoils along this
  ////  line to figure out how wide (and where) the scintillator should be.
  ////	   
  
  double pr = sqrt(SQR(newRecoilMomentum.vec.x) + SQR(newRecoilMomentum.vec.y));
  double pz = newRecoilMomentum.vec.z;
  
  // compute slope in r-z plane
  // tv_ prefix denotes 2-component vec (we use the first two components
  // of a Three vector to hold r,z)
  //  i.e.
  // (x,y,z) = (r,z,<dummy>)
 	    
  double pmt_angle_rads = pmt_angle_degrees * 2.0 * PI / 360.0;
  
  Three tv_targetCenter(0.0,
			target_length / 2.0, 
			0.0 /* dummy */);
  
  Three tv_aux1(target_to_pmt_dist * sin(pmt_angle_rads),
		target_to_pmt_dist * cos(pmt_angle_rads), 
		
		0.0 /* dummy */);
  
  Three tv_point1 = tv_targetCenter + tv_aux1;
  // point where line from target center intersects PMT plane at right angle
  
  Three tv_aux2(0.0,target_to_pmt_dist / cos(pmt_angle_rads), 0.0 /* dummy */);
  
  Three tv_point2 = tv_targetCenter + tv_aux2;
  // point where PMT plane intersects z-axis
  
  Three tv_PMTplane = tv_point2 - tv_point1;
  
  double PMTplaneSlope = tv_PMTplane.vec.x / tv_PMTplane.vec.y;
  // actually the r component over the z component
  
  double PMTrIntercept = (-PMTplaneSlope) * target_to_pmt_dist / cos(pmt_angle_rads);
  
  // TEST
  
  // 	    	tv_point1.Print();
  //	tv_point2.Print();
  //	printf("m, b = %f, %f \n",PMTplaneSlope,PMTrIntercept);
  
  // END TEST
  
  
  // now find point where recoil intersects plane
  // See page 45 of my CEBAF notes
  // PCF 27 JAN 96
  
  
  Three tv_bendPoint( sqrt(SQR(bendPoint.vec.x) + SQR(bendPoint.vec.y)),
		      bendPoint.vec.z,
		      0.0);
  
  
  double u = 
    (PMTrIntercept + (PMTplaneSlope*(tv_bendPoint.vec.y)) - (tv_bendPoint.vec.x))
    / (pr - (PMTplaneSlope*pz));
  
  
  Three tv_recoilP( pr,pz,0.0);
  
  Three tv_PMTIntersectionPoint = tv_bendPoint + tv_recoilP*u;
  
  
  //tv_PMTIntersectionPoint.Print();
  
  // compute distance from point where PMT plane intersects (z) axis
  
  
  double PMTdist = (tv_PMTIntersectionPoint - tv_point2).Mag();


  ////
  //// Make plots.
  //// Also, a couple of last minute calculations for the plots.
  ////
  
  //  Generated events	
  
  //  let's play with stuff; MORE HISTOGRAMS!
  z=tv_PMTIntersectionPoint.vec.y;
  rho=tv_PMTIntersectionPoint.vec.x;
  theta=(atan(rho/z)*180)/PI;
  /*  bf1(52,theta,1.0); */
  // this is the generated version
  //  CzorS 03JUN96
  int recoiled;
  if(theta>THETA_MIN && theta<THETA_MAX && currentRecoilMomentum>minRecoil_momentum)
  {
    recoiled=1;
    //    printf("made recoil cut\t%f\t%f\n",theta,currentRecoilMomentum);
  }
  else
  {
    recoiled=0;
    //    printf("didn't make recoil cut\t%f\t%f\n",theta,currentRecoilMomentum);
  }
  
  reject_flag=glass_sim(vertex,recoiled);
  
  bf2(101,p_recoil,angle,1.0); // momentum/angle
  bf1(106,angle,1.0);          // plot recoil angle
  
  bf1(102,target.mag3(),1.0);      // total fermi momentum
  bf1(103,target.vec.space.x,1.0); // x,y,z comps of fermi
  bf1(104,target.vec.space.y,1.0);
  bf1(105,target.vec.space.z,1.0);
  
  /*  bf1(10,beam.vec.space.x,1.0); // beam momenta
  bf1(11,beam.vec.space.y,1.0);
  bf1(12,beam.vec.space.z,1.0);
  bf1(13,beamTheta,1.0);
  
  bf1(107,vertex.vec.x,1.0);
  bf1(108,vertex.vec.y,1.0);
  bf1(109,vertex.vec.z,1.0);
  */
  if(!(isnand(recoil.vec.t-currentRecoilEnergy)))
    bf1(110,recoil.vec.t - currentRecoilEnergy,1.0); // energy loss plot
  
  float newRecoilCosTheta = (newRecoilMomentum * zAxis) / newRecoilMomentum.Mag(); 
  double newRecoilAngle = acos(newRecoilCosTheta) * 360.0 / 6.28318; 
  
  bf1(111,newRecoilAngle,1.0);
  if(!(isnand(currentRecoilMomentum))&&!(isnand(newRecoilAngle)))
    bf2(112,currentRecoilMomentum,newRecoilAngle,1.0); 
  // energy angle with m.s. and e. loss
  bf1(113,ms_theta,1.0);                   
  // net scatter angle in radians
  
  double recoilKE = recoil.vec.t - recoilMass;
  double recoilKEafterLoss = currentRecoilEnergy - recoilMass;
  
  /*  bf1(150,recoilKE,1.0);
  bf1(151,recoilKEafterLoss,1.0);
  bf2(152,recoilKE,angle,1.0);
  bf2(153,recoilKEafterLoss,newRecoilAngle,1.0);
  */
  if(gEscapedTarget)
    {
      //      bf1(154,PMTdist,1.0);
      //bf2(155,tv_PMTIntersectionPoint.vec.x,tv_PMTIntersectionPoint.vec.y,1.0);  
      // r,z of recoil in PMT plane
    };
  
  ////
  //// See if this event was accepted or not, if it's okay so far.
  ////   	
  
  rejected = DoCebafAcceptance(currentRecoilMomentum);  	

  double cos_theta_recon;
  
  if(reject_flag==0)     //changed 30SEP96 CzorS
    {
      // Accepted events

      b_r_CM=into_rest_frame_of(reconstructed_CM,smeared_product);
      b_r_bat=into_rest_frame_of(smeared_g_in_glass[config-1],smeared_product);
      
      cos_theta_recon=costheta_between(b_r_CM,b_r_bat);
      bf1(12012,cos_theta_recon,1.0);
      bf2(12310,cos_theta_recon,cos(qtheta_rad),1.0);

      bf2(12300,true_beam_E,smeared_product_E,1.0);

      if(!isnand(smeared_product_mass))
	bf1(12202,smeared_product_mass,1.0);

      bf2(10101,p_recoil,angle,1.0); // momentum/angle
      bf1(10106,angle,1.0);          // plot recoil angle
		  
      bf1(10102,target.mag3(),1.0);      // total fermi momentum
      bf1(10103,target.vec.space.x,1.0); // x,y,z comps of fermi
      bf1(10104,target.vec.space.y,1.0);
      bf1(10105,target.vec.space.z,1.0);
      
      /*      bf1(10010,beam.vec.space.x,1.0); // beam momenta
      bf1(10011,beam.vec.space.y,1.0);
      bf1(10012,beam.vec.space.z,1.0);
      bf1(10013,beamTheta,1.0);
      */
      /*      bf1(10107,vertex.vec.x,1.0);
      bf1(10108,vertex.vec.y,1.0);
      bf1(10109,vertex.vec.z,1.0);
      */
      bf1(10110,recoil.vec.t - currentRecoilEnergy,1.0); // energy loss plot 
      bf1(10111,newRecoilAngle,1.0);
      bf2(10112,currentRecoilMomentum,newRecoilAngle,1.0);// energy angle with m.s. and e. loss
      bf1(10113,ms_theta,1.0);  // net scatter angle in radians
      /*      bf1(10150,recoilKE,1.0);
      bf1(10151,recoilKEafterLoss,1.0);
      bf2(10152,recoilKE,angle,1.0);
      bf2(10153,recoilKEafterLoss,newRecoilAngle,1.0);
      */
      bf1(10087,pt_squared,1.0);
      bf1(12001,minus_t_recoil,1.0);
      bf1(12003,minus_t_gammas,1.0);

      bf1(12011,cos(qtheta_rad),1.0);

      bf2(12101,(-(t_reac_true)),pt_squared,1.0);

      if(beta>0.001)  //arbitrary; this is really slow, though
	bf1(10089,(.66667/beta),1.0); //time to go 20cm


      if(gEscapedTarget)
	{
	  //	  bf1(10154,PMTdist,1.0);
	  //bf2(10155,tv_PMTIntersectionPoint.vec.x,tv_PMTIntersectionPoint.vec.y,1.0);  
	  // r,z of recoil in PMT plane
	};
      
      //  let's play with stuff; MORE HISTOGRAMS!
      rho=tv_PMTIntersectionPoint.vec.x;
      z=tv_PMTIntersectionPoint.vec.y;
      theta=(atan(rho/z)*180)/PI;
      if(gEscapedTarget)
	//	bf1(10052,theta,1.0);
      //  CzorS
      g_nAccepted++;
    }
  else // the event didn't make it
    {
      // process rejects here
      
    };	
  
  g_nTotal++;

 
  void *dummy_pointer=NULL;
#define GROUP_GEO_MC_TRACKS 100
  dummy_pointer=data_addGroup(ih,BUFF_SIZE,
     GROUP_GEO_MC_TRACKS,0,
     (unsigned long)(sizeof(geo_tracks_rec_t)+((NUM_TRACKS-1)*sizeof(geo_track_t))) );
  
  memcpy(dummy_pointer,reaction_tracks,
     (sizeof(geo_tracks_rec_t)+((reaction_tracks->ntracks-1)*sizeof(geo_track_t))) );
  
  //int data_write(int fd,const void* event);   /*automatically regenerates the CRC word*/

  
  if(buffer!=NULL)
    data_write(fileno(buffer),ih);
  
} // end of doStuff

////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

Three FindVertex(Three inEntryPoint, Three inMomentum)
{
   ////// determine the vertex in a cylindrical target.
   ////// Inputs are where the beam intersects the upstream
   ////// side of the target, and the beam momentum.
   //////
   //////  Output is the vertex position vector.  It
   //////  is equally likely to lie anywhere on the line
   //////  that connects the target entry point and
   //////  the target exit point.
   //////
   //////  PCF 21 JAN 96
   
   
   Three exitPoint = FindExitPoint(inEntryPoint, inMomentum);
      
   double s = (exitPoint.vec.z - inEntryPoint.vec.z) / inMomentum.vec.z;
   
   // Pick a place uniformly along this line; that's the vertex
   
   double x = s * unfRand();  // random number between 0 and s
   
   Three vertex = inEntryPoint + inMomentum * x;
   return vertex;
   
}


///////////////////////////////////////////////////////////////////////

Three FindExitPoint(Three inStartPoint, Three inMomentum)
{

        ////////////////// this routine, given a starting 3-vector
        ////////////////// and momentum, will determine the point
        ////////////////// of exit from the target cylinder.
        //////////////////
        ////////////////// geom params are in "targetRegion.h"
        //////////////////
        ////////////////// assumes (0,0,0) is center of upstream
        //////////////////  face of the target
        //////////////////
        //////////////////  PCF 19 JAN 96
        
        
   
   Three currentPosition,outExitPoint;
   float t = 2.0;
 
    
   // first find a point outside of the target
   
   int done = 0;
   
   do
   {
     currentPosition = inStartPoint + inMomentum*t; 
     if ( PointInTarget(&currentPosition) )
     { 
       done = 0;
       t *= 2.0;
     }
     else
     {
       done = 1;
     }
    } while (!done);
    
   // now do binary search to find exit point
   
   done = 0;
   
   float tin = 0;
   float tout = t;
   
   int loops=0;
   do
   {
     /* this loop is hanging the program.  I'm going to fix it.
	Craig Steffen 03FEB99 */
     t = tin + (tout - tin) / 2.0;
     currentPosition = inStartPoint + inMomentum*t;
    
     if(PointInTarget(&currentPosition))
     {
       tin = t;
     }
     else
     {
       tout = t;
     };
     
     if ( ((tout - tin) / tout) < 0.00001 ) done = 1;
   
     loops++;
   } while(!done&&loops<10000);
   
   outExitPoint = inStartPoint + inMomentum*tout;
   
   return outExitPoint;
   
}

int PointInTarget(Three *inPointToTest)
{
	// returns 1 if point is in target
	// 0 if not
	//
	// assumes center of upstream face of target is (0,0,0)
	
	double r_coord = sqrt( SQR(inPointToTest->vec.x) + SQR(inPointToTest->vec.y) );
	double z_coord = inPointToTest->vec.z;
	
	if ((r_coord <= target_radius) && (z_coord >= 0) && (z_coord <= target_length))
	  return 1;
	  
	return 0;
}

///////////////////////////////////////////////////////////////////////

void printInfo(void)
{

  int qa;
  float ratio;
  float pct_Acc = 100.0 * (float)g_nAccepted / (float)g_nTotal;
  float pct_MT =  100.0 * (float)g_nMissingTarget / (float)(gParameters.codes.nIterations);
  float pct_RHE = 100.0 * (float)g_nRecoilEnHi / (float)g_nTotal;
  float pct_RLE = 100.0 * (float)g_nRecoilEnLow / (float)g_nTotal;
  
  printf("\n\n");
  printf("CEBAF exited successfully\n\n");
  printf("%d events were generated.\n",gParameters.codes.nIterations);
  printf("%5.1f %% of these missed the target.\n\n",pct_MT);
  printf("%d total events interacted in target.\n\n",g_nTotal);
  //  printf("%5.1f %% had recoil energy too high.\n", pct_RHE);
  printf("%5.1f %% had recoil energy too low.\n\n", pct_RLE);
  //  printf("%5.1f %% pass all cuts.\n", pct_Acc);
  for(qa=MAX_CUTS-1;qa>0;qa=qa-1)
    {
      if(reject_reason_totals[qa]>0)
	{
	  ratio=(((float)reject_reason_totals[qa]*100)/((float)(g_nTotal)));
	  printf("%d (%4.3f %s) events were rejected\n",reject_reason_totals[qa],ratio,"%");
	  printf("\tbecause %s.\n",reasons[qa]);
	}
    }
  printf("\n\nEvents making all Cebaf cuts:\t\t%d\t%4.3f%s\n",total_accepted,
	 ((float)(total_accepted*100)/(float)(g_nTotal))," %");
  printf("\nEvents making partial acceptance:\t%d\t%4.3f%s\n"
	 ,partial_cuts_accepted,
	 ((float)(partial_cuts_accepted*100)/(float)(g_nTotal))," %");
  
  
  printf("\n\n");


  //tabulation file output:
  if(tabulator!=NULL)
    {
      fprintf(tabulator,"//////////////////////////////////////////////////\n\n");
      //    printf("\n%s\n",&beam_type[0]);



      fprintf(tabulator,"%s + %s  -->  %s + %s\n\n",beam_type,target_type,product_type,recoil_type);
      fprintf(tabulator,"Run number  %d\n",run_number);
      fprintf(tabulator,"Electron beam energy : %3.1f GeV\n",gElectronBeamE);

      
      fprintf(tabulator,"%d events were generated.\n",gParameters.codes.nIterations);
      fprintf(tabulator,"%d total events interacted in target.\n\n",g_nTotal);
      fprintf(tabulator,"%5.1f %% had recoil energy too low.\n\n", pct_RLE);  
      for(qa=MAX_CUTS-1;qa>0;qa=qa-1)
	{
	  if(reject_reason_totals[qa]>0)
	    {
	      ratio=(((float)reject_reason_totals[qa]*100)/((float)(g_nTotal)));
	      fprintf(tabulator,"%d (%4.3f %s) events were rejected\n",reject_reason_totals[qa],ratio,"%");
	      fprintf(tabulator,"\tbecause %s.\n",reasons[qa]);
	    }
	}
      fprintf(tabulator,"\n\nEvents making all Cebaf cuts:\t\t%d\t%4.3f%s\n",total_accepted,
	      ((float)(total_accepted*100)/(float)(g_nTotal))," %");
      fprintf(tabulator,"\nEvents making partial acceptance:\t%d\t%4.3f%s\n"
	      ,partial_cuts_accepted,
	      ((float)(partial_cuts_accepted*100)/(float)(g_nTotal))," %");
      
      
      fprintf(tabulator,"\n\n");
      fprintf(tabulator,"/////////////////////////////////////////////////\n");
    }
}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////

void SaveThemHistos(void)
{
  bbsave(0,"cebaf.bbook");
}
////////////////////////////////////////////////////////////////////////

int DoCebafAcceptance(double theRecoilMomentum)
{
  int accept = 0; // 0 good, 1 bad
  
  if (theRecoilMomentum < minRecoil_momentum)
    { 
      accept = 1;
      g_nRecoilEnLow++;
    };
  
  if (theRecoilMomentum > maxRecoil_momentum)
    {
      accept = 1;
      g_nRecoilEnHi++;
    };
  
  
  return accept;
}

////////////////////////////////////////////////////////////////////////

void PrintWelcomeMessage(void)
{
  printf("\n\nCEBAF E-94-016 target region monte carlo\n");
  printf("========================================\n");
  printf("\n\n");
  
  printf("Electron beam energy : %3.1f GeV\n",gElectronBeamE);
  printf("Max beam angle       : %3.1f theta C's\n",gParameters.beam.bremThetaMax);
  printf("Target radius        : %4.2f cm\n",target_radius);
  printf("Target length        : %4.2f cm\n",target_length);
  printf("Recoil momentum cut  : %4.2f GeV/c\n",minRecoil_momentum);
  
  printf("========================================\n");
  printf("\n\nGenerating...\n\n");
}

////////////////////////////////////////////////////////////////////////

int glass_sim(Three vertx,int recoil_passed)
{
  int qindex,inglass,rej[MAX_CUTS],least_sep=0,checker=0,dflag,dflag2,tempcount;
  int rej_this_time[MAX_CUTS];
  float X[MAX_GAMMA],Y[MAX_GAMMA],rho[MAX_GAMMA];
  float E_gamma[MAX_GAMMA],separation[MAX_GAMMA];
  float rho2;
  Three vertices[MAX_GAMMA],mom[MAX_GAMMA];
  lgdGeom_t impact_block;
  vector3_t reac_vertex;
  double distance;
  for(qindex=0;qindex<MAX_CUTS;qindex++)
    {
      rej[qindex]=0;
      rej_this_time[qindex]=0;
    }
  //config
  /*  if(config==2)
    tempcount=2;
  else
    if(config==3)
      tempcount=3; */
  for(qindex=0;qindex<config;qindex++)
    {
      /*          printf("%f  %f  %f %f ",g_in_glass[qindex].vec.space.x,g_in_glass[qindex].vec.space.y,g_in_glass[qindex].vec.space.z,g_in_glass[qindex].vec.t); */

      distance=(GLASS_TO_TARGET_Z-vertx.vec.z)/g_in_glass[qindex].vec.space.z; 
      /*      distance=135.0; */
      /*      fprintf(stdout,"qindex=%d\n",qindex); */
            mom[qindex]=Three(g_in_glass[qindex].vec.space); 
	    /*      mom[qindex].vec.x=g_in_glass[qindex].vec.space.x;
      mom[qindex].vec.y=g_in_glass[qindex].vec.space.y;
      mom[qindex].vec.z=g_in_glass[qindex].vec.space.z;
      fprintf(stdout,"qindex=%d\n",qindex); */
      vertices[qindex]=(distance*mom[qindex])+vertx;
      /*      printf("%f  %f  %f\n",vertices[qindex].vec.x,vertices[qindex].vec.y,vertices[qindex].vec.z); */
      E_gamma[qindex]=g_in_glass[qindex].vec.t;
      bf1(10067,E_gamma[qindex],1.0);
      if(E_gamma[qindex]<GAMMA_E_MIN)
	{
	  rej[5]=1;
	  if(rej_this_time[5]==0)
	    {
	      reject_reason_totals[5]++;
	      rej_this_time[5]=1;
	    } 
	}
      X[qindex]=vertices[qindex].vec.x;
      Y[qindex]=vertices[qindex].vec.y;
      reac_vertex.x=(float)(vertices[qindex].vec.x);
      reac_vertex.y=(float)(vertices[qindex].vec.y);
      reac_vertex.z=(float)(vertices[qindex].vec.z);
      /* Craig Steffen's Kludge */
      impact_block.row=(int)((reac_vertex.x+56.0)/4.0);
      impact_block.col=(int)((reac_vertex.y+56.0)/4.0);
      rho2=sqrt(SQR(reac_vertex.x)+SQR(reac_vertex.y));
      impact_block.space.x=(float)((impact_block.row*4)+2);
      impact_block.space.y=(float)((impact_block.col*4)+2);
      impact_block.space.z=0.0;
      impact_block.blockSize.x=4.0;
      impact_block.blockSize.y=4.0;
      impact_block.blockSize.z=45.0;
      
      if((rho2>5.66)&&(rho2<51.0)&&(g_in_glass[qindex].vec.space.z>0.0)){
	/* considered accepted for 620 element detector */
       	impact_block.channel=((impact_block.row*28)+impact_block.col);
	if(impact_block.channel>783||impact_block.channel<0){
	  printf("found channel %d, for row %d column %d; exiting.\n",impact_block.channel,impact_block.row,impact_block.col);
	}
	impact_block.hole=0;
      }
      else{
	bf1(12004,g_in_glass[qindex].vec.t,1.0);
	  bf2(12005,acos(g_in_glass[qindex].vec.space.z/g_in_glass[qindex].vec.t),g_in_glass[qindex].vec.t,1.0);
	impact_block.channel=-1;
	impact_block.row=-1;
	impact_block.col=-1;
	impact_block.hole=1;
      }
      
      
      /* inglass=lgdGlobalPositionToGeomType(reac_vertex,&impact_block); */
      /* inglass=lgdMPSPositionToGeomType(reac_vertex,&impact_block);*/
      /* printf("x=%f, y=%f, z=%f, r=%d, c=%d, hole=%d.\n",impact_block.space.x,impact_block.space.y,impact_block.space.z,impact_block.row,impact_block.col,impact_block.hole); */
      
      
      if(impact_block.hole==LGDGEOM_TRUE)
	{
	  rej[1]=1;
	  if(rej_this_time[1]==0)
	    {
	      reject_reason_totals[1]++;
	      rej_this_time[1]=1;
	    }     
	}
      else
	{
	  if(lgdGeomIsAnEdge(impact_block.channel)==LGDGEOM_TRUE)
	    {
	      rej[6]=1;
	      if(rej_this_time[6]==0)
		{
		  reject_reason_totals[6]++;
		  rej_this_time[6]=1;
		}
	    }
	}
      rho[qindex]=sqrt((X[qindex]*X[qindex])+(Y[qindex]*Y[qindex]));
      bf2(10058,X[qindex],Y[qindex],1.0);
      bf1(10059,rho[qindex],1.0);
    }
  /*  printf("\n"); */
  if(recoil_passed==0)
    {
      rej[2]=1;
      if(rej_this_time[2]==0)
	{
	  reject_reason_totals[2]++;
	  rej_this_time[2]=1;
	}
    }

  /* this is where the separation logic was. */

  if(separation[least_sep]<MIN_GAMMA_SEP)
    {    
      rej[4]=1;
      if(rej_this_time[4]==0)
	{
	  reject_reason_totals[4]++;
	  rej_this_time[4]=1;
	}
    }
  /*  rej[0]+=rej[1];
  for(qindex=3;qindex<MAX_CUTS;qindex++)
    rej[0]+=rej[qindex]; */
  for(qindex=1;qindex<MAX_CUTS;qindex++)
    rej[0]+=rej[qindex]; 
  
  
  //////////////////////////////////////////////////////////
  //           reject table, 15JUL96, CzorS
  // rej[1]<>0 when >0 gammas missed the glass
  // rej[2]<>0 when recoil was not detected
  // rej[3]<>0 
  // rej[4]<>0 when two photns are closer than 8 cm
  // rej[5]<>0 when one or more gammas failed energy cuts
  // rej[6]<>0 when gamma landed in an edge block
  // 
  // rej[0] is the sum of the above, =0 when event accepted
  //
  //////////////////////////////////////////////////////////
  
  
  
  //gamma radius, raw
  {
    //    bf1(10073,rho[dflag],1.0);
    //  bf1(10074,rho[dflag2],1.0);
  }
  
  //gamma radius, partial cuts, defined here
  if(rej[1]==0&&rej[4]==0&&rej[5]==0&&rej[6]==0)
    {
      bf1(10071,rho[dflag],1.0);
      bf1(10072,rho[dflag2],1.0);
      partial_cuts_accepted++;
      bf1(10094,(-(t_reac_true)),1.0);      
    }

  if(rej[0]==0)
    {
      bf1(10085,(-(t_reac_true)),1.0);
      bf1(10081,rho[dflag],1.0);
      bf1(10082,rho[dflag2],1.0);
      total_accepted++;
      //config
      if(glassfile!=NULL)
	fprintf(glassfile,"%d\n",config);
      for(qindex=0;qindex<config;qindex++)
	{
	  //	printf("%d\t%f\t%f\n",qindex,X[qindex],Y[qindex]);
	  bf2(10068,X[qindex],Y[qindex],1.0);
	  bf1(10069,rho[qindex],1.0);
	  bf1(10080,separation[least_sep],1.0);
	  
	  if(glassfile!=NULL)
	    {
	      fprintf(glassfile,"%f\t%f\t%f\t%f\n",g_in_glass[qindex].vec.t,g_in_glass[qindex].vec.space.x,g_in_glass[qindex].vec.space.y,g_in_glass[qindex].vec.space.z);
	      fprintf(glassfile,"%f\t%f\t%f\n",vertices[qindex].vec.x,vertices[qindex].vec.y,vertices[qindex].vec.z);
	    }
	}
      if(glassfile!=NULL)
	fprintf(glassfile,"\n");
    }
  return(rej[0]);
}




void printHelp(void)
{
  fprintf(stderr, "Usage:  Cebaf -i<Chees file> -r<run #> -c<config #> -t<tab file>\n");
  fprintf(stderr, "\n<Chees file> is a Chees configuration file name with the extension .xrdb.\n");
  fprintf(stderr, "<run #> is the run number.  At this point, this is used to choose the\n\t size of the lgd.\n");
  fprintf(stderr, "<config #>The number of photons in the final state.  This should be 2 or 3.\n");
  fprintf(stderr, "<tab file> is a file name of the file where you are tabulating results.\n\tThis is on append mode, so use the same file more than once.\n");
  fprintf(stderr, "  -h\t\tPrint this message\n\n");
  fprintf(stderr, "An example command line is:\n");
  fprintf(stderr, "Cebaf -iphi.xrdb -r1000 -c3 -ttabulate.out\n\n");
}

  
main(int argc, char *argv[])
{

  char *inputFilename = NULL;
  char *flag = NULL;
  char *run_char;
  char *outputFilename =NULL;
  char *tabulateFilename=NULL;
  
  reasons[1]="one or more gammas missed the glass";
  reasons[2]="recoil was not detected";
  reasons[4]="two photons were too close together";
  reasons[5]="one or more gammas failed energy cuts";
  reasons[6]="one or more gammas landed in edge blocks";



  /* glassfile=fopen("glassvecs.txt","w"); */
  /*  if(!glassfile)
    {
      printf("problem opening glassvecs.txt");
      exit(29);
    }*/
  // First some command line BS
  // Not very interesting, skip down a bit.
  //
  // Get input filename from command line.
  
  int qz;
  for(qz=0;qz<MAX_CUTS;qz++)
    reject_reason_totals[qz]=0;
  
  for(int iarg=1; iarg < argc; iarg++)
    if(*argv[iarg] == '-')
      {
        flag = argv[iarg]+1;
	switch(*flag)
          {
	  case 'i':
	    inputFilename = flag+1;
	    break;
	  case 'o':
	    outputFilename= flag+1;
	    //	    outputFilename= argv[++iarg];
	    if(!(buffer=fopen(outputFilename,"w"))){
	      fprintf(stderr,"can't open output file... exiting\n");
	      exit(1);
	    }
	    break;
	  case 't':
	    tabulateFilename=flag+1;
	    if(!(tabulator=fopen(tabulateFilename,"a")))
	      {
		printf("tabulator file cannot be opened\n");
		exit(1);
	      }
	    break;
	  case 'c':
	    config=atoi(flag+1);
	    break;
	  case 'r':
	    run_char=flag+1;
	    run_number=atoi(run_char);
	    break;
	  case 'h':
	    fprintf(stderr, "\n");
	    printHelp();
	    exit(0);
	  default:
	    fprintf(stderr, "\nInvalid command line specification\n");
	    printHelp();
	    exit(EXIT_FAILURE);
	  }
      }

  if(!(inputFilename))  /* bastards didn't give us an input file */ 
    {
      printf("\n\nPlease specify an input file, e.g.,\n\n");
      printHelp();
      exit(EXIT_FAILURE);
    }
  if(run_number==0)
    {
      printf("run number is 0, this is invalid\n\n");
      exit(18);
    }
  printf("Run number is \t%d\n\n",run_number);

  ih=(itape_header_t*)malloc(BUFF_SIZE);
  if(!ih)
    {
      fprintf(stderr,"Failed to initialize buffer, exiting.\n");
      exit(2);
    }

  CutChees(inputFilename, &gParameters, &gParticleDBPtr);  // init CHEES
  
  gElectronBeamE = gParameters.beam.p;
  
  PrintWelcomeMessage();
  InitThemHistos();
  switch(gParameters.codes.geometry){
  case GEOM_CEBAF_328:
    lgdGeomSetup(-200);
    break;
  case GEOM_CEBAF_620:
    lgdGeomSetup(-100);
    break;
  default:
    fprintf(stderr,"Geomtry must be CEBAF_328 or CEBAF_620\n");
    exit(2);
  }
  
    // GENERATE EVENTS:
  
  for (int event = 0; event < gParameters.codes.nIterations; event++)
    {
        if(!(event % 500 ))
      {
	printf("%d\r",event); 
	//printf("%d\r\n",event); 
	fflush(stdout);
      }
      
      doStuff();    // the heart of the program
      
    };
  
  ///
  /// Done genreating all events, clean up.
  ///
  
  SaveThemHistos();   
  
  ScrubChees(); // clean up CHEES
  
  
  /* --------------------------------------------------------- */

  printInfo();

  fclose(buffer);
  fclose(tabulator);
  fclose(glassfile);

  return 0;
}


////////////////////////////////////////////////////////////////////////

/* This is a dummy comment */

Four smear_4v(Four invec)
{
  Four retvec;
  double xslope, yslope, zslope;
  xslope = invec.vec.space.x / invec.vec.t;
  yslope = invec.vec.space.y / invec.vec.t;
  zslope = invec.vec.space.z / invec.vec.t;

  double sigma = 
    ENERGY_SMEAR_FLOOR  + (ENERGY_SMEAR_STAT /(sqrt(invec.vec.t)));

  sigma *= invec.vec.t;

  short done = 0;
  while (!done)   /* pick a new energy, making sure it's greater than 0 */
  {

    double newEnergy =  
     GetGaussianDist(invec.vec.t, sigma);

    if (newEnergy > 0.0)
    {
     retvec.vec.t = newEnergy;
     done = 1;
    };

  } /* end while loop */

  retvec.vec.space.x = xslope * retvec.vec.t;
  retvec.vec.space.y = yslope * retvec.vec.t;
  retvec.vec.space.z = zslope * retvec.vec.t;  

  retvec=smearPosition_4v(retvec);

  return retvec;
  
}

Four smearPosition_4v(Four invec){
  double x,y,xslope,yslope;

  x=(invec.vec.space.x/invec.vec.space.z)*(GLASS_TO_TARGET_Z+10.0);
  y=(invec.vec.space.y/invec.vec.space.z)*(GLASS_TO_TARGET_Z+10.0);

  x = GetGaussianDist(x, LGD_X_SMEAR_VALUE);
  y = GetGaussianDist(y, LGD_Y_SMEAR_VALUE);
  xslope = (x)/(GLASS_TO_TARGET_Z+10.0);
  yslope = (y)/(GLASS_TO_TARGET_Z+10.0);
  invec.vec.space.z = invec.vec.t / sqrt(SQR(xslope) + SQR(yslope) + 1);
  invec.vec.space.x = xslope * invec.vec.space.z;
  invec.vec.space.y = yslope * invec.vec.space.z;
  return invec;
}
