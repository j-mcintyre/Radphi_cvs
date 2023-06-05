/*
 * mcwrap - adapted by R.T.Jones for Radphi
 *          from cwrap by S. Teige
 *
 * This program reads an input file that specifies a series of
 * two-body or specific three-body decays for products of the
 * photoproduction reaction.  Currently three-body decays with
 * non-uniform Dalitz plots for omega, eta and etaprime are
 * supported.  The output is in the Radphi itape format ready
 * for input to the Gradphi GEANT simulator.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>


/* Misc Radphi */

#include <itypes.h>
#include <disData.h>
#include <particleType.h>

/* Radphi I/O */

#include <dataIO.h>
#include <disIO.h>
#include <tapeIO.h>
#include <tapeData.h>

#include <mctypes.h>

/* control variables */

static int flags[100];
static float emax=8.0;
static float emin=6.0;

/* Itape building, I/O, interfaces etc. */

void Usage(void);
void setupItape(itape_header_t **event);
void makeMCraw(itape_header_t *event, int istat);
Particle_t giveType( int id );
int giveCharge(int id );
int getone_(int *iflag, float *, float * , char *);
void mcgive_(float *, float *, float *, int *id, int *iflag, int *itot);


/* cernlib randoms generator */

void rluxgo_(int *lux, int *seed, const int *k1, const int *k2);
int lux_level = 3;
int lux_k1 = 0;
int lux_k2 = 0;

/* End of function declarations */

/* I/O variables */

#define BUFSIZE 300000
FILE *fpOutRaw=NULL;
int fileNo=0;

/* Misc event variables */

itape_header_t *event=NULL;
int runNo=9999;
int spillNo=1;
int eventNo=0;

void Usage(void)
{
  fprintf(stderr,"Usage :\n");
  fprintf(stderr,"cwrap <Options>\n");
  fprintf(stderr,"Options are :\n");
  fprintf(stderr,"\t-h\t\tPrint this message.\n");
  fprintf(stderr,"\t-D\t\tDebug output (default: no debug output)\n\n");

  fprintf(stderr,"\t-r#\t\tUse # as random number seed.\n");
  fprintf(stderr,"\t-l#\t\tUse # as minimum tagged photon energy. (default: 6 GeV)\n");
  fprintf(stderr,"\t-u#\t\tUse # as maximum tagged photon energy. (default: 8 GeV)\n\n");

  fprintf(stderr,"\t-c<filename>\tOutput to <filename>.itape in c mode.(default: no output)\n");
  fprintf(stderr,"\t-a\t\tWrite intermediate state particles to output (c output mode only,\n");
  fprintf(stderr,"\t\t\tdefault: write only stable particles)\n\n");
  fprintf(stderr,"\t-F<filename>\tOutput to <filename>.mcdat in FORTRAN mode.(default:, no output)\n");
  fprintf(stderr,"\t-q\t\tQuiet mode. No generation messages.\n");
}

int main(int argc, char *argv[])

{
  int iarg;
  char *argptr;
  int index;
  int randSeed;  
  char *ftnFile=NULL;
  char *outputFile=NULL;
  int istat;
  char fname[1024];
  char ftnname[1024];

  randSeed = getpid() + time(NULL);
  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'u':
	emax = atof(argptr+1);
	break;
      case 'l':
	emin = atof(argptr+1);
	break;
      case 'r':
	randSeed = atoi(argptr+1);
	break;
      case 'q':
        flags[0]=1;
	break;
      case 'o':
	break;
      case 'h':
	Usage();
	return;
      case 'F':
	ftnFile = argptr+1;
	flags[1]=1;
	break;
      case 'c':
	outputFile = argptr+1;
	flags[2]=1;
	break;
      case 'D':
        flags[3]=1;
	break;
      case 'a':
        flags[4]=1;
	break;
      default:
	fprintf(stderr,"Unknown argument : %s\n");
	Usage();
	return;
      }
    }
  if(!flags[0]) fprintf(stderr,"Random Seed : %d\n",randSeed);
  if(!flags[0]) fprintf(stderr,"Photon energies from %f to %f\n",emin,emax);

  //srandom(randSeed); -- I don't know what this is, but it does not work -- rtj 3-3-2020
  rluxgo_(&lux_level, &randSeed, &lux_k1, &lux_k2);

  /* 
   * Start the output up 
   */
  if(flags[1]){
    sprintf(ftnname,"%s.mcdat",ftnFile);
    if(!flags[0]) fprintf(stderr,"Writing FORTRAN format raw MC to %s\n",ftnname);
  }
  if(flags[2]){
    sprintf(fname,"%s.itape",outputFile);
    if(!flags[0]) fprintf(stderr,"Writing c format raw MC to %s\n",fname);
    fpOutRaw = fopen(fname,"w");
  }
  /* 
   *Loop over requested number of events 
   */
  istat=2;
  while(istat != 1){
    istat=getone_(flags,&emin,&emax,ftnname);
    if(istat != 1){
      eventNo++;
      /*
       * Print something out now and then so I know its alive
       */
      if(!(eventNo%1000) && (eventNo!=0)){
	fprintf(stderr,"%d \r",eventNo);
	fflush(stdout);
      }
      if(flags[2]){
	setupItape(&event);
	/*
	 * Write all events.
	 */
	makeMCraw(event,istat);
	if(fpOutRaw){
	  data_clean(event);
	  data_write(fileno(fpOutRaw),event);
	}
      }
    }
  }  
  /*
   * Cleanup the output files
   */
  if(flags[2]){
    if(fpOutRaw){
      data_flush(fileno(fpOutRaw));
      fclose(fpOutRaw);
    }
  }
  return 0;
}

void setupItape(itape_header_t **event)
{
  if(!(*event))
    *event = (itape_header_t*)malloc(BUFSIZE);
  
  data_newItape(*event);

  (*event)->runNo = runNo;
  (*event)->spillNo = eventNo/100+1;
  (*event)->eventNo = eventNo;
  (*event)->trigger = 1;

  /* Set the event type eventually */

}

void makeMCraw(itape_header_t *event, int istat)

{
  mc_event_t *mc_event=NULL;
  mc_vertex_t *mc_vertex=NULL;
  int index;
  static int initialSize, finalSize;
  static int icount;
  int itot;
  float beamdump[4];
  float pl[5][50];
  float xyzv[3][50];
  int id[50];
  int iflag[50];

  if(istat == 2){
    /*
     * have a new type of event must recalculate its size
     */
    mcgive_(beamdump,&pl[0][0],&xyzv[0][0],id,iflag,&itot); 
    icount=0;
    for(index=0;index<49;index++){
      if(flags[4]==1){
	/*
	 * have requested intermediate state particles calculate size appropriately.
	 */
	if((iflag[index]!=10000)&&(iflag[index]!=0)){
	  icount++;
	}
      }
      else{
	/*
	 * Are writing only stable final state particles.
	 */
	if((iflag[index]==50000)||(iflag[index]==70000)){
	  icount++;
	}
      }
    }
    initialSize = sizeof_mc_event_t(2);
    finalSize = sizeof_mc_vertex_t(icount);
  }

  /* 
   * initialize counters in the MC structure 
   */
  mc_event = data_addGroup(event,BUFSIZE,GROUP_MC_EVENT,0,initialSize);
  mc_vertex = data_addGroup(event,BUFSIZE,GROUP_MC_EVENT,1,finalSize);
  /*
   * get Monte-Carlo information
   */
  mcgive_(beamdump,&pl[0][0],&xyzv[0][0],id,iflag,&itot);
  /*
   * Add content to the MC structure 
   */
  mc_event->type=0;
  mc_event->weight=1;
  mc_event->nfinalVert=0;
  mc_event->npart=0;
  /* 
   * primary vertex
   */
  mc_vertex->kind = 0;
  mc_vertex->origin.t = 0;
  mc_vertex->origin.space.x = xyzv[0][0];
  mc_vertex->origin.space.y = xyzv[1][0];
  mc_vertex->origin.space.z = xyzv[2][0];
  mc_vertex->npart = 0;
  mc_event->nfinalVert++;
  /* 
   * beam vector
   */
  mc_event->part[0].momentum.space.x = beamdump[0];
  mc_event->part[0].momentum.space.y = beamdump[1];
  mc_event->part[0].momentum.space.z = beamdump[2];
  mc_event->part[0].momentum.t = beamdump[3];
  mc_event->part[0].kind = 1;
  mc_event->npart++;
  /*
   * final state particles
   */
  icount=0;
  for(index=0;index<itot;index++){
    /*
     * Look for particles to go to the output list.
     * First decide what output list we are making.
     */
    if(flags[4]!=0){
      /*
       * Are writing intermediate and final state particles.
       */
      if(iflag[index]!=10000){
	/*
	 * Write everything except the beam.
	 */
	mc_vertex->npart++;
	mc_vertex->part[icount].momentum.space.x=pl[0][index];
	mc_vertex->part[icount].momentum.space.y=pl[1][index];
	mc_vertex->part[icount].momentum.space.z=pl[2][index];
	mc_vertex->part[icount].momentum.t=pl[3][index];
	mc_vertex->part[icount].kind=giveType(id[index]);
	icount++;
      }
    }
    else{
      /*
       * Are writing only stable particles,
       */
      if((iflag[index]==50000)||(iflag[index]==70000)){
	/*
	 * have found a stable particle.
	 */
	mc_vertex->npart++;
	mc_vertex->part[icount].momentum.space.x=pl[0][index];
	mc_vertex->part[icount].momentum.space.y=pl[1][index];
	mc_vertex->part[icount].momentum.space.z=pl[2][index];
	mc_vertex->part[icount].momentum.t=pl[3][index];
	mc_vertex->part[icount].kind=giveType(id[index]);
	icount++;
      }
    }
  }
}
Particle_t giveType(int id)
{
  /*
   * Translation from fffnmc particle types to GEANT particle types
   */
  if(id==8){
    return(Neutron);
  }
  if(id==9){
    return(Proton);
  }
  if(id==11){
    return(KLong);
  }
  if(id==12){
    return(KShort);
  }
  if(id==13){
    return(KMinus);
  }
  if(id==14){
    return(KPlus);
  }
  if(id==15){
    return(PiPlus);
  }
  if(id==16){
    return(PiMinus);
  }
  if(id==17){
    return(Pi0);
  }
  if(id==18){
    return(Electron);
  }
  if(id==19){
    return(Positron);
  }
  if(id==20){
    return(Neutrino);
  }
  if(id==21){
    return(Gamma);
  }
  if(id==22){
    return(MuonPlus);
  }
  if(id==23){
    return(MuonMinus);
  }
  if(id==25){
    return(Eta);
  }
  if(id==26){
    return(Rho0);
  }
  if(id==27){
    return(omega);
  }
  if(id==28){
    return(EtaPrime);
  }
  if(id==31){
    return(PhiMeson);
  }
  return(0);
}
int giveCharge(int id)
{
  /*
   * Translation fffnmc id to charge
   */
  if(id==8){
    return(0);
  }
  if(id==9){
    return(1);
  }
  if(id==13){
    return(-1);
  }
  if(id==14){
    return(1);
  }
  if(id==15){
    return(1);
  }
  if(id==16){
    return(-1);
  }
  if(id==21){
    return(0);
  }
  return(0);
}
