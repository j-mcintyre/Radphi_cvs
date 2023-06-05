/*
 * lgdtune - stand-alone program to tune gain coefficients on LGD adcs
 *           so as to optimize the experimental resolution of some peak
 *           in the invariant mass spectra.
 *
 * - rev 1.2 - Richard Jones, May 4 2001.
 *
 * Modified to facilitate the production of calibration on large data sets.
 *    1. Introduced the option -d<dropfile>.  The dropfile contains the
 *	 internal state of the calibrator at the end of the last iteration
 *	 pass.  A sequence of dropfiles can be loaded into a subsequent job
 *	 at initialization time by use of the -a<dropfile> option.  Use of
 *	 the -d option automatically forces -i1 (single-iteration job).
 *    2. Introduced the option -a<dropfile>.  The dropfile must have been
 *	 created by a previous lgdtune job that carried the -d<dropfile>
 *	 option.  By using the -a option, lgdtune picks up where the former
 *	 job left off and can continue to accumulate calibration statistics
 *	 and solve for tune coefficients.  A check is made that the dropfile
 *	 is compatible with the current job, that is, that it was produced
 *	 using the same set of gain coefficients as is being used by the
 *	 present job.  If the job starts from a non-default tune (using the
 *	 -f option) then the -f<file> must appear somewhere among the options.
 *	 One can use multiple -a options, in which case the sequence of
 *	 dropfiles is accumulated as if the former processing had taken place
 *	 in a single lgdtune job.  Each dropfile encountered must be compatible
 *	 with this job, in the sense discussed above.  Use of the -a option
 *	 automatically forces -i1 (single-iteration job).
 *
 * - rev 1.1 - Richard Jones, April 26 2001.
 *
 * Modified the selection of events and mesons to be included in calibration
 *    1. Suppress events with more than 3 clusters to control combinatoric
 *	 background under the peak being optimized.
 *    2. Make a cut on the geometric mean energy of two photons that make up
 *	 a pi0.  This is because large values for E1 E2 correspond to small
 *	 opening angles between the photons.  The minimum two-cluster distance
 *	 imposed by cluster finding and cleanup produces a large bias in the
 *	 pi0 mass twards high values for pairs whose minimum opening angle
 *	 puts the average cluster pair distance within this radius.  It is
 *	 controlled by the following macro.
 */
#define PAIR_MIN_OPENING_ANGLE 0.09
#define MAX_CLUSTER_MULTIPLICITY 3
/*
 *
 * - original code - Richard Jones, March 2 2001.
 *
 * Differences from lgdCalibrate:
 *    The code is based on S. Teige's lgdCalibrate.  In particular the
 *    user interface was kept the same.   The major differences are these.
 *    1. The chi-square is now defined as the difference in mass-squared
 *       between measured and physical value.  In lgdCalibrate it was defined
 *       between measured and fitted cluster energy, which required a fit of
 *       cluster sets to kinematic hypotheses.  I wanted to avoid the fit.
 *    2. A Lagrange multiplier has been added to the chi-square in order to
 *       remove the bias in the corrected peak position.  Without this term
 *       the fit always places the peak some fraction of a sigma below the
 *       true mass.  This is correct behavior mathematically, but is not what
 *       we want.  This effect is described in detail in the next section.
 *    3. All cluster finding, corrections, pruning, etc are done in standard
 *       library routines.  This means that the mass spectra that come out
 *       from lgdtune will be the same ones as are observed by subsequent
 *       analyses that use the tuned gain coefficients.
 *    4. Matrix math is done by cernlib routines, instead of private code.
 *       The linear algebra package in cernlib is fairly extensive, well
 *       debugged and has plenty of flexibility for variable precision.
 *    5. The output data file contains new information: the final gain
 *       coefficients, the correction factors from the overall run, and
 *       the correction factors from the last step, to show the degree
 *       of convergence.
 *    6. Comments.
 *
 * Method:
 *    The basic idea is that well-tuned gain coefficients should minimize
 *    the observed width of prominent peaks in invariant mass spectra.
 *    Candidate mesons are pi0, eta, omega, and maybe etaprime; the only
 *    requirement is that their physical widths be small compared to the
 *    experimental resolution.  Starting from the existing calibration from
 *    the database, the algorithm calculates correction factors to the gain
 *    of each tube which approximately minimize the observed peak width.
 *    After several passes through the calibration loop, the correction
 *    factors should converge to 1, indicating that the optimum has been
 *    found.  The optimization function is a constrained chi-squared:
 *
 *          F{gc_i} = Sum_n [ mm^2 - m_n^2 ]^2 +
 *             lambda Sum_n [ mm^2 - m_n^2 ]			(1)
 *
 *    where mm is the true mass of the meson and m_n is the mass seen
 *    in event n (or occurrence n allowing for more than one per event) 
 *    which is a function of the gain parameters gc_i.  The first term
 *    (quadratic in mass-squared) is the ordinary chi-squared function
 *    in mass-squared.  I have not worried about rescaling it by putting
 *    a sigma in the denominator because it just affects F by an overall
 *    scale.  The second term in (1) is there to enforce the constraint
 *    that the peak be centered on the correct mass after the fit.
 *    NOTE:  This does NOT follow automatically from the first term!
 *           The reason for this is that, there are actually two ways for
 *           the fit to reduce the width of the peak:
 *           (a) adjust the gains to make it narrower (good)
 *           (b) leave the relative gains fixed and just reduce the
 *               overall scale factor toward zero (bad)
 *    By doing (b) the fit reduces both the mean and the sigma, the mean
 *    costing and the sigma paying the first term.  They reach equilibrium
 *    somewhere, depending on the experimental width of the peak.  But
 *    if you use just the quadratic term alone, you will get a systematic
 *    downward bias in the peak of the corrected spectra that you will have
 *    to take out somehow.  The second term uses a Lagrange multiplier
 *    lambda to constrain the bias to zero.
 *
 * Details of algorithm:
 *    The mass function m_n is written as a function of the gc_i as follows:
 *
 *          m_n^2 = Sum_p [ Sum_ij [ 2 gc_i gc_j E_i E_j
 *                                  * (1-cos(theta_12)) ] ]	(2)
 *
 *    where the sum is over all cluster-pairs p and over all blocks i in
 *    cluster 1 and blocks j of cluster 2 in the pair.  The angle theta_12
 *    between the two clusters is considered to be independent of the
 *    tube gains, an approximation.  By this formula we can see that the
 *    standard mass-squared is obtained when the gc_i=gc_j=1.  One also can
 *    see that the gradient of m_n^2 in the space of the {gc_i} is just the
 *    standard mass-squared of the pair that block participates in times the
 *    energy fraction E_i/E that the given block contributes to its cluster.
 *
 *          partial_i(m_n^2) = (m_n_p)^2 * f_i			(3)
 *
 *    where m_n_p is what appears inside the Sum_p in equation (2).
 *    Formally the solution for gain correction factors gc_i is given by
 *
 *          gc_i  = 1 + Cinverse D - lambda Cinverse L         	(4)
 *
 *    where the matrix Cinverse is the inverse of the curvature matrix,
 *
 *          C_ij = Sum_n [ partial_i(m_n^2) partial_j(m_n^2) ]	(5)
 *          D_i  = Sum_n [ partial_i(m_n^2) ( mm^2 - m_n^2 )	(6)
 *          L_i  = Sum_n [ partial_i(m_n^2) ]			(7)
 *
 *    At the end one has to solve for lambda to satisfy the zero-bias
 *    constraint.
 *
 *    One has to handle the matrix Cinverse carefully because it is subject
 *    to parasitic effects that spoil the convergence of the solution.  These
 *    effects come about because Cmatrix is subject to statistical errors
 *    which make the measured eigenvalues vary somewhat from the true ones.
 *    Suppose because of these errors that one of the eigenvalues of Cmatrix
 *    comes out very small or zero.  In such a large matrix this is almost
 *    unavoidable, but they don't make much difference to Cmatrix.  Remember
 *    that the real symmetric Cmatrix can be written in terms of its eigen-
 *    values c_m and eigenvectors Vm as
 *
 *                C_ij  = sum_m ( c_m Vm_i Vm_j )       (2)
 *
 *    So small c_m do not affect C_ij very much.  But consider its inverse
 *    which can be written in terms of the same eigenvalues and eigenvectors.
 *
 *          Cinverse_ij = sum_m ( (1/c_m) Vm_i Vm_j)    (3)
 *
 *    Now one small eigenvalue c_m can dominate the entire structure of
 *    Cinverse.  The region of the detector most vulnerable to this sort
 *    of instability is the outer periphery where the statistics on an
 *    individual block are poor and fluctuations are the largest.  Even
 *    for very large statistics, if enough iterations are run one will
 *    eventually find a small eigenvalue and convergence will be lost.
 *    The only effective way to control it is to impose a cutoff on the
 *    smallest eigenvalue c_m allowed in the expansion (2).  That is why
 *    I work with an eigenvector decomposition of Cmatrix and Cinverse
 *    instead of directly calculating the matrix inverse.
 */

#define MIN_CMATRIX_EIGENVALUE 10
#define HBOOK 8000000
#define NDIM 784


#include <stdio.h>
#include <stdlib.h>

#include <disData.h>
#include <dataIO.h>
#include <itypes.h>
#include <eventType.h>
#include <triggerType.h>
#include <lgdCluster.h>
#include <lgdGeom.h>
#include <lgdUtil.h>
#include <particleType.h>
#include <makeHits.h>
#include <tapeData.h>
#include <umap.h>
#include <math.h>
#include <cernlib.h>

int main(int argc, char *argv[]);

/* LAPACK driver for finding eigenvals/vects for real symmetric matrices */

void dspev_(const char *jobz,  /* 'N' for eigenvalues only, 'V' for both    */
            const char *uplo,  /* 'U' for upper triangular A, 'L' for lower */
            const int *n,      /* order of matrix A                         */
            double *A,         /* input real symmetric matrix, overwritten  */
            double *w,         /* output eigenvalues in ascending order     */
            double *z,         /* output eigenvectors in order of w         */
            const int *ldz,    /* leading dimension of z                    */
            double *work,      /* work space, 3*n elements                  */
            int *info);        /* exit status, 0 on success                 */

/* internal state of the calibrator */

static double Cmatrix[NDIM][NDIM];
static double Dvector[NDIM];
static double Lvector[NDIM];
static double massbias;
static float gain[NDIM];
static float gc[NDIM];

/* spectral decomposition tables for Cmatrix */

double eigenval[NDIM];
double eigenmode[NDIM][NDIM];
int channels, dof;

/* parameter defaults for this run */

int nChannels=NDIM;
int nIterations=1;
float energyCut=4.0;
int mesonCut=7;
float windowCut=0.02;
float base_gc=1.0;
char *gcFile=NULL;
char *dropFile=NULL;
char *preFile[9999];
int npreFiles=0;

/* Need to make cc output files... */

static float cc_out[NDIM];
static char *ccOut=NULL;

struct hbk{float ia[HBOOK];} pawc_;


void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  lgdCalibrate <options> file1 file2 ...\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-i<n>\t\titerate tuning loop <n> times,");
  fprintf(stderr," default %d\n", nIterations);
  fprintf(stderr,"\t-m<n>\t\toptimize on meson type <m>,");
  fprintf(stderr," default %d (%s)\n", mesonCut,ParticleType(mesonCut));
  fprintf(stderr,"\t-w<GeV>\t\tset mass window width to <GeV>,");
  fprintf(stderr," default %f\n", windowCut);
  fprintf(stderr,"\t-c<fact>\tpreset coefficients to <fact>,");
  fprintf(stderr," default %f\n", base_gc);
  fprintf(stderr,"\t-E<GeV>\t\trequire total energy <GeV> in LGD,");
  fprintf(stderr," default %f\n", energyCut);
  fprintf(stderr,"\t-f<file>\ttake initial gains from <file>,");
  fprintf(stderr," default none\n");
  fprintf(stderr,"\t-M<max>\t\tanalyze maximum of <max> events,");
  fprintf(stderr," default none\n");
  fprintf(stderr,"\t-d<file>\tsave final internal state to <file>,");
  fprintf(stderr," default none\n");
  fprintf(stderr,"\t-a<file>\taccumulate internal state from <file>,");
  fprintf(stderr," default none\n");
  fprintf(stderr,"\t-o\tsave cc info into cc.out\n");
  fprintf(stderr,"\t-O<file>\tsave cc info into <file>\n");
  fprintf(stderr,"\t-h\t\tprint this message\n");
}

meson_search_list_t *shopping_list(void)
{
  meson_search_list_t *myList=malloc(sizeof_meson_search_list_t(99));
  int stage=0;

/* Define the meson search algorithm for this run here.
 * Any meson in the list can be used as a tune target.
 * Only a few mesons make sense to use as tune targets,
 * (pi0, eta, omega, maybe phi, etaprime) because their
 * natural widths should be much smaller than the experimental
 * mass resolution.  The tune works best if the peak in the
 * mass plot is not sitting on top of a lot of background.
 */

  myList->try[stage].meson = Pi0;
  myList->try[stage].daughter[0] = Gamma;
  myList->try[stage].daughter[1] = Gamma;
  myList->try[stage].width = 0.100;
  ++stage;
     
  myList->try[stage].meson = Eta;
  myList->try[stage].daughter[0] = Gamma;
  myList->try[stage].daughter[1] = Gamma;
  myList->try[stage].width = 0.150;
  ++stage;

  myList->try[stage].meson = omega;
  myList->try[stage].daughter[0] = Pi0;
  myList->try[stage].daughter[1] = Gamma;
  myList->try[stage].width = 0.200;
  ++stage;

  myList->try[stage].meson = PhiMeson;
  myList->try[stage].daughter[0] = Eta;
  myList->try[stage].daughter[1] = Gamma;
  myList->try[stage].width = 0.200;
  ++stage;

  myList->try[stage].meson = EtaPrime;
  myList->try[stage].daughter[0] = Gamma;
  myList->try[stage].daughter[1] = Gamma;
  myList->try[stage].width = 0.200;
  ++stage;

  myList->try[stage].meson = f0_980;
  myList->try[stage].daughter[0] = Pi0;
  myList->try[stage].daughter[1] = Pi0;
  myList->try[stage].width = 0.300;
  ++stage;

  myList->try[stage].meson = a0_980;
  myList->try[stage].daughter[0] = Eta;
  myList->try[stage].daughter[1] = Pi0;
  myList->try[stage].width = 0.300;
  ++stage;

  myList->stages = stage;
  return myList;
}

void gcSolver()
{
  const int ldim=NDIM;
  static double A[NDIM*(NDIM+1)/2];
  double *work;
  int info;
  int channel[NDIM];
  int i,j,k;

/* sequeeze out the null columns */

  channels=0;
  for(i=0;i<nChannels;i++){
    if(Cmatrix[i][i]>0){
      channel[channels++]=i;
    }
  }

/* pack Cmatrix into upper-diagonal form */

  k=0;
  for(j=0;j<channels;j++){
    int jb=channel[j];
    for(i=0;i<=j;i++){
      int ib=channel[i];
      A[k++]=Cmatrix[ib][jb];
    }
  }

/* use LAPACK routine to solve eigenvalue problem */

  work=malloc(sizeof(double)*3*channels);
  dspev_("V","U",&channels,A,
          eigenval,(double *)eigenmode,&ldim,work,&info);
  free(work);
  if(info<0){
    fprintf(stderr,"dspev error: invalid argument in position %d\n",-info);
    exit(1);
  }
  else if(info>0){
    fprintf(stderr,"dspev error: convergence error %d\n",info);
    exit(1);
  }
  else {

/* Instead of computing Cinverse directly and doing matrix products, we
 * chose to do an orthogonal transform on all vectors into the frame in
 * which C is diagonal.  Express Dvector and Lvector in eigenmode expansion.
 */
    double Dmode[NDIM];
    double Lmode[NDIM];
    double LCinvL,LCinvD;
    double lambda;
    int mode,block,iblock;
    dof=0;
    for(mode=0;mode<channels;mode++){
      Dmode[mode]=0;
      Lmode[mode]=0;
      for(iblock=0;iblock<channels;iblock++){
        block=channel[iblock];
        Dmode[mode]+=eigenmode[mode][iblock]*Dvector[block];
        Lmode[mode]+=eigenmode[mode][iblock]*Lvector[block];
      }
      if(eigenval[mode]>MIN_CMATRIX_EIGENVALUE){
        ++dof;
      }
    }

/* Solve for lambda = (L Cinverse D  -  massbias)/( L Cinverse L) */

    LCinvL=0;
    LCinvD=0;
    for(mode=0;mode<channels;mode++){
      if(eigenval[mode]>MIN_CMATRIX_EIGENVALUE){
        LCinvL+=Lmode[mode]*Lmode[mode]/eigenval[mode];
        LCinvD+=Lmode[mode]*Dmode[mode]/eigenval[mode];
      }
    }
    lambda=(LCinvD-massbias)/LCinvL;

/* Load the results into the gc vector */

    for(block=0;block<nChannels;block++){
      gc[block]=1;
    }
    for(iblock=0;iblock<channels;iblock++){
      block=channel[iblock];
      for(mode=0;mode<channels;mode++){
        if(eigenval[mode]>MIN_CMATRIX_EIGENVALUE){
          int block=channel[iblock];
          float strength=Dmode[mode]-lambda*Lmode[mode];
          gc[block]+=strength*eigenmode[mode][iblock]/eigenval[mode];
        }
      }
    }
  }
}

void setupHbook(void)
{
  int index;
  char str[1024];

  hlimit(HBOOK);

  hbook1(1,"total energy",200,0,8,0);
  hbook1(2,"nClusters",10,-0.5,9.5,0);
  hbook2(3,"hits",28,-0.5,27.5,28,-0.5,27.5,0);  
  hbook2(4,"hits in clusters",28,-0.5,27.5,28,-0.5,27.5,0);  
  hbook2(5,"gain tunes",28,-0.5,27.5,28,-0.5,27.5,0);
  hbook1(6,"last gc",100,0,2,0);
  hbook1(7,"Cmatrix eigenvalues",100,0,100,0);
  hbook1(8,"two-gamma mass spectrum",300,0,1.2,0);
  for(index=2;index<9;index++){
    hbook1(8+index,"two-gamma mass spectrum",300,0,1.2,0);
  }
  hbook2(9,"gain tunes, adjusted channels",28,-0.5,27.5,28,-0.5,27.5,0);
  for(index=1;index<99;index++){
    char *pname=ParticleType(index);
    if (strcmp(pname,"unknown") == 0) {
      sprintf(str,"%s spectrum",pname);
      hbook1(1000+index,str,100,0,2*ParticleMass(index),0);
    }
  }
  hbook1(25,"gain tunes",100,0,2,0);

  /* cc ouput file histos */
  hbook1(50,"cc ouput values",250,0,0.005,0);
  hbook1(51,"delta cc",250,-0.0005,0.0005,0);

  hbook1(100,"total energy, all events",200,0,8,0);
}

int configureIt(itape_header_t *event)
{
  int actualRun;
  int index;
  int index2;
  int jpoint;
  float x;

  meson_search_list_t *myList;

  if(setup_makeHits(event)){
    fprintf(stderr,"Problem with setup_makeHits; exiting");
    return;
  }
  myList=shopping_list();
  if(setup_makeMesons(event,myList)){
    fprintf(stderr,"Problem with setup_makeMesons; exiting");
    return;
  }
  free(myList);
  return(0);
}

void processEvent(itape_header_t *event)
{

  float Etot;
  float Etot_clus;
  lgd_hits_t *lgdHits;
  lgd_clusters_t *clusters;
  lgd_hits_t *clusterHits;
  lgd_hits_t *scaledHits;
  photons_t *photons;
  mesons_t *mesons;
  int size;
  int row,col;
  int index;

  if((event->eventType==EV_DATA)&&(event->trigger==TRIG_DATA)){
    if(unpackEvent(event,BUFSIZE)){
      fprintf(stderr,"processEvent: unpacking error!\n");
      return;
    }
    if(make_lgd_hits(event,NULL,640)){
      fprintf(stderr,"processEvent: make_lgd_hits error!\n");
      return;
    }
    lgdHits = data_getGroup(event,GROUP_LGD_HITS,0);
    if((lgdHits==NULL)||(lgdHits->nhits==0))
      return;

    Etot=0;
    for(index=0;index<lgdHits->nhits;index++){
      int channel=lgdHits->hit[index].channel;
      /* Apply gain adjustment here */
      lgdHits->hit[index].energy*=gain[channel];
      Etot += lgdHits->hit[index].energy;
      row=channel/28;
      col=channel-(row*28);
      hf2(3,col,row,lgdHits->hit[index].energy);
    }

    hf1(100,Etot,1);

    if(Etot<energyCut){
      return;
    }
    hf1(1,Etot,1);

    if(make_lgd_clusters(event,NULL,100)){
      fprintf(stderr,"processEvent: make_lgd_clusters error!\n");
      return;
    }
    clusters = data_getGroup(event,GROUP_LGD_CLUSTERS,0);
    clusterHits = data_getGroup(event,GROUP_LGD_CLUSTER_HITS,0);
    if((clusters->nClusters==0)||(clusterHits->nhits==0))
      return;

    Etot_clus=0;
    for(index=0;index<clusterHits->nhits;index++){
      int channel=clusterHits->hit[index].channel;
      Etot_clus += clusterHits->hit[index].energy;
      row=channel/28;
      col=channel-(row*28);
      hf2(4,col,row,clusterHits->hit[index].energy);
    }

    if(lgd_cluster_cleanup(event)){
      return;
    }
    hf1(2,clusters->nClusters,1);

    if(make_photons(event,NULL,100)){
      fprintf(stderr,"processEvent: make_photons error!\n");
      return;
    }

/* To calculate derivatives below we will need to have the hit
 * energies scaled as fractions of their cluster energy.
 */
    size=sizeof_lgd_hits_t(clusterHits->nhits);
    scaledHits=memcpy(malloc(size),clusterHits,size);
    for(index=0;index<clusters->nClusters;index++){
      float energy=clusters->cluster[index].energy;
      int iblock=clusters->cluster[index].firstBlock;
      int i;
      for(i=0;i<clusters->cluster[index].nBlocks;i++){
        scaledHits->hit[iblock++].energy/=energy;
      }
    }

    photons = data_getGroup(event,GROUP_PHOTONS,0);
    for(index=0;index<clusters->nClusters;index++){
      int index2;
      for(index2=index+1;index2<clusters->nClusters;index2++){
        float pair_mass2;
        vector3_t p;
        float p_0;
        p_0=photons->photon[index].energy+
            photons->photon[index2].energy;
        p.x=photons->photon[index].momentum.x+
            photons->photon[index2].momentum.x;
        p.y=photons->photon[index].momentum.y+
            photons->photon[index2].momentum.y;
        p.z=photons->photon[index].momentum.z+
            photons->photon[index2].momentum.z;
        pair_mass2=p_0*p_0-p.x*p.x-p.y*p.y-p.z*p.z;
        hf1(8,sqrt(pair_mass2),1);

	if (clusters->nClusters<9){
         hf1(8+clusters->nClusters,sqrt(pair_mass2),1);
	}

      }
    }

    if(make_mesons(event,NULL,100)){
      fprintf(stderr,"processEvent: make_mesons error!\n");
      free(scaledHits);
      return;
    }
    mesons = data_getGroup(event,GROUP_MESONS,0);

/* Do only events with MAX_CLUSTER_MULTIPLICITY events or fewer
 */ if(photons->nPhotons > MAX_CLUSTER_MULTIPLICITY){
      free(scaledHits);
      return;
    }
    for(index=0;index<mesons->nMesons;index++){
      Particle_t type=mesons->meson[index].type;
      float mass=mesons->meson[index].mass;
      hf1(1000+type,mass,1);
      if(type==mesonCut){
        float mmass=ParticleMass(mesonCut);
        float mass_diff=mmass-mass;
        if(fabs(mass_diff)<windowCut){

/* Find the list of all gammas that belong to this mother meson
 */       int glist[10],ngammas=0;
          int stack[999],pstack=0;
          stack[pstack++]=index;
          while(pstack>0){
            int this=stack[--pstack];
            int d;
            for(d=0;d<mesons->meson[this].ndecay;d++){
              if(mesons->meson[this].dtype[d]==Gamma)
                glist[ngammas++]=mesons->meson[this].daughter[d];
              else
                stack[pstack++]=mesons->meson[this].daughter[d];
            }
          }

/* For a two-photon decay, require a certain minimum opening angle
 */       if(ngammas==2){
            int g1=glist[0];
            int g2=glist[1];
            float gmeanE=sqrt(photons->photon[g1].energy
                             *photons->photon[g2].energy);
            float gmeanEmax=mmass/PAIR_MIN_OPENING_ANGLE;
            if(gmeanE>gmeanEmax){
              ngammas=0;
            }
          }

/* Find the gradient of the mass-squared wrt gain of each block present
 */       if(ngammas>1){
            double grad[NDIM];
            int block[NDIM];
            int blocks;
            int iblock;
            int ig1,ig2;
            float mmass2=mmass*mmass;

/* The overall mass-squared is just the sum of mass-squared of each pair.
 * Sum the mass-squared for all pairs that each photon participates in.
 */         float mass2[10];
            memset(mass2,0,sizeof(float)*10);
            for(ig1=0;ig1<ngammas;ig1++){
              int g1=glist[ig1];
              for(ig2=ig1+1;ig2<ngammas;ig2++){
                int g2=glist[ig2];
                float pair_mass2;
                vector3_t p;
                float p_0;
                p_0=photons->photon[g1].energy+
                    photons->photon[g2].energy;
                p.x=photons->photon[g1].momentum.x+
                    photons->photon[g2].momentum.x;
                p.y=photons->photon[g1].momentum.y+
                    photons->photon[g2].momentum.y;
                p.z=photons->photon[g1].momentum.z+
                    photons->photon[g2].momentum.z;
                pair_mass2=p_0*p_0-p.x*p.x-p.y*p.y-p.z*p.z;
                mass2[ig1]+=pair_mass2;
                mass2[ig2]+=pair_mass2;
              }
            }

/* The gradient of the mass-squared function for each block is just the
 * above mass-squared sum over all pairs involving the given cluster
 * times the fraction of that cluster's energy provided by that block.
 */         blocks=0;
            for(ig1=0;ig1<ngammas;ig1++){
              int g1=glist[ig1];
              int c1=photons->photon[g1].cluster;
              int i;
              iblock=clusters->cluster[c1].firstBlock;
              for(i=0;i<clusters->cluster[c1].nBlocks;i++){
                block[blocks]=scaledHits->hit[iblock].channel;
                grad[blocks]=scaledHits->hit[iblock].energy*mass2[ig1];
                ++blocks;
                ++iblock;
              }
            }

/* Accumulate the results in Dvector and Cmatrix
 */         for(iblock=0;iblock<blocks;iblock++){
              int b1=block[iblock];
              int iblock2;
              for(iblock2=iblock;iblock2<blocks;iblock2++){
                int b2=block[iblock2];
                int bgt=(b1>b2)?b1:b2;
                int blt=(b1>b2)?b2:b1;
                Cmatrix[blt][bgt]+=grad[iblock]*grad[iblock2]/(mmass2*mmass2);
              }
              Dvector[b1]+=grad[iblock]*(1-mass*mass/mmass2)/mmass2;
              Lvector[b1]+=grad[iblock]/mmass2;
            }
            massbias+=(1-mass*mass/mmass2);
          }
        }
      }
    }
    free(scaledHits);
  }
}

int doAnEvent(int fd)
{
  static itape_header_t *event=NULL;
  int ret;
  tapeHeader_t *hdr=NULL;
  static int currentRun=0;

  if(!event)
    event = malloc(BUFSIZE);
  
  ret=evt_data_read(fd,event,BUFSIZE);
  switch(ret){      
  case DATAIO_OK:
    switch(event->type){
    case TYPE_TAPEHEADER:
      hdr=(void *)event;
      return(doAnEvent(fd));		
    case TYPE_ITAPE:
      if(event->runNo!=currentRun){
	if(configureIt(event)){
	  fprintf(stderr,"Configuration failed!\n");
	  return(1);
	}
	currentRun = event->runNo;
      }
      processEvent(event);
      return(0);
    default:
      /* Oh no, no idea what this is.... */
      fprintf(stderr,"doAnEvent: got event of unknown type : %d\n",event->type);
      return(doAnEvent(fd));
    }
  case DATAIO_EOF:
    /* Tape marker, more data to come */
    fprintf(stderr,"doAnEvent: EOF on tape\n");
    return(doAnEvent(fd));
  case DATAIO_EOT:
    /* Real EOF tape or file */
    return(1);
  case DATAIO_ERROR:
    fprintf(stderr,"doAnEvent: data_read() returned DATAIO_ERROR,");
    fprintf(stderr," Skipping event\n");
    return(1);
  case DATAIO_BADCRC:
    fprintf(stderr,"doAnEvent: data_read() returned DATAIO_BADCRC,");
    fprintf(stderr," Skipping event\n");
    return(1);
  }
  return(1);
}

int main(int argc, char *argv[])
{
  int iarg;
  char *argptr;
  int fdIn;
  int maxEvents=0;
  int nEvents;
  int index;
  int iteration;
  char defoutFile[]="lgdtune.out";
  char *outFile=defoutFile;
  int Use_ccOut=0;
  char defccOut[]="cc.out";
  char *ccOut=defccOut;
  int row=0,col=0;

  for(iarg=1;iarg<argc;iarg++){
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'm':
	mesonCut = atoi(++argptr);
	break;
      case 'w':
	windowCut = atof(++argptr);
	break;
      case 'c':
	base_gc = atof(++argptr);
	break;
      case 'f':
	gcFile = (char *)++argptr;
	break;
      case 'E':
	energyCut = atof(++argptr);
	break;
      case 'M':
	maxEvents = atoi(++argptr);
	break;
      case 'i':
	nIterations = atoi(++argptr);
	break;
      case 'd':
	dropFile = (char *)++argptr;
	break;
      case 'a':
	preFile[npreFiles++] = (char *)++argptr;
	break;
      case 'o':
	Use_ccOut=1;
	break;
      case 'O':
	Use_ccOut=1;
	ccOut = (char *)++argptr;
      case 'h':
	Usage();
	exit(1);
      default:
	fprintf(stderr,"Unknown argument: %s\n",argv[iarg]);
	Usage();
	exit(1);
      }
    }
  }
  
  setupHbook();

  if(gcFile){
    FILE *fpin=fopen(gcFile,"r");
    int col,row,chan;
    if(fpin==NULL){
      fprintf(stderr,"lgdtune error: input file %s",gcFile);
      fprintf(stderr," not found\n");
      exit(1);
    }
    for(index=0;index<nChannels;index++){
      fscanf(fpin,"%d %d %d %f %f",
                   &chan,&col,&row,&gain[index],&gc[index]);
      if(chan!=index){
        fprintf(stderr,"lgdtune error: channel %d out of order",chan);
        fprintf(stderr," in input file %s\n",gcFile);
        exit(1);
      }
    }
    fclose(fpin);
  }
  else {
    for(index=0;index<NDIM;index++){
      gain[index] = base_gc;
    }
  }

  if (nIterations != 1) {
    if (npreFiles != 0) {
      fprintf(stderr,"lgdtune error: -i option incompatible with -a\n");
      exit(1);
    }
    if (dropFile != 0) {
      fprintf(stderr,"lgdtune error: -i option incompatible with -d\n");
      exit(1);
    }
  }

/* Here we fill the Cmatrix, Dvector, ... from data
 * which (hopefully) are rich in the designated meson
 */
  for(iteration=0;iteration<nIterations;iteration++){
    int ipreFile;
    memset(Cmatrix,0,sizeof(double)*NDIM*NDIM);
    memset(Dvector,0,sizeof(double)*NDIM);
    memset(Lvector,0,sizeof(double)*NDIM);
    massbias=0;
    hreset(0," ");
    nEvents=0;

/* read in any initialization from dropfiles generated by previous jobs */

    for(ipreFile=0;ipreFile<npreFiles;ipreFile++){
      static double lCmatrix[NDIM][NDIM];
      double lDvector[NDIM];
      double lLvector[NDIM];
      double lmassbias;
      float lgain[NDIM];
      float lgc[NDIM];
      int nev;
      FILE *fpin=fopen(preFile[ipreFile],"r");
      if(fpin==NULL){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," not found\n");
        exit(1);
      }
      if(fread(lCmatrix,NDIM*sizeof(double),NDIM,fpin)<NDIM){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," returned error during read\n");
        exit(1);
      }
      if(fread(lDvector,sizeof(double),NDIM,fpin)<NDIM){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," returned error during read\n");
        exit(1);
      }
      if(fread(lLvector,sizeof(double),NDIM,fpin)<NDIM){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," returned error during read\n");
        exit(1);
      }
      if(fread(&lmassbias,sizeof(double),1,fpin)<1){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," returned error during read\n");
        exit(1);
      }
      if(fread(lgain,sizeof(float),NDIM,fpin)<NDIM){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," returned error during read\n");
        exit(1);
      }
      if(fread(lgc,sizeof(float),NDIM,fpin)<NDIM){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," returned error during read\n");
        exit(1);
      }
      if(fread(&nev,sizeof(int),1,fpin)<1){
        fprintf(stderr,"lgdtune error: input file %s",preFile[ipreFile]);
        fprintf(stderr," returned error during read\n");
        exit(1);
      }
      for(index=0;index<NDIM;index++){
        int ii;
        if(lgain[index]!=gain[index]){
          fprintf(stderr,"lgdtune error: gains from %s",preFile[ipreFile]);
          fprintf(stderr," are incompatible with those for this job\n");
          exit(1);
        }
        for(ii=0;ii<NDIM;ii++){
          Cmatrix[index][ii]+=lCmatrix[index][ii];
        }
        Dvector[index]+=lDvector[index];
        Lvector[index]+=lLvector[index];
      }
      massbias+=lmassbias;
      nEvents+=nev;
      fclose(fpin);
    }

    for(iarg=1;iarg<argc;iarg++){
      if(*argv[iarg] != '-'){
	if( (fdIn = evt_open(argv[iarg],"r")) == 0){
	  fprintf(stderr,"Failed to open %s for reading\n",argv[iarg]);
	  exit(1);
	}
	else{
	  fprintf(stderr,"Reading from %s\n",argv[iarg]);
	  while( (maxEvents ? maxEvents > nEvents : 1) && !doAnEvent(fdIn)){
            int i;
            nEvents++;
            for(i=10000;i<nEvents;i*=10){}
            i/=10;
            if(nEvents/i*i == nEvents)
	      fprintf(stderr,"%d\n",nEvents);
	  }
	  if(fdIn){
	    evt_close(fdIn);
	    fdIn = 0;
	  }
	}      
      }
    }

/* generate a dropfile if requested */

    if(dropFile!=NULL){
      FILE *fpout=fopen(dropFile,"w");
      if(fpout==NULL){
        fprintf(stderr,"lgdtune error: drop file %s",preFile[ipreFile]);
        fprintf(stderr," could not be opened for writing\n");
        exit(1);
      }
      if(fwrite(Cmatrix,NDIM*sizeof(double),NDIM,fpout)<NDIM){
        fprintf(stderr,"lgdtune error: drop file %s",dropFile);
        fprintf(stderr," returned error during write\n");
        exit(1);
      }
      if(fwrite(Dvector,sizeof(double),NDIM,fpout)<NDIM){
        fprintf(stderr,"lgdtune error: drop file %s",dropFile);
        fprintf(stderr," returned error during write\n");
        exit(1);
      }
      if(fwrite(Lvector,sizeof(double),NDIM,fpout)<NDIM){
        fprintf(stderr,"lgdtune error: drop file %s",dropFile);
        fprintf(stderr," returned error during write\n");
        exit(1);
      }
      if(fwrite(&massbias,sizeof(double),1,fpout)<1){
        fprintf(stderr,"lgdtune error: drop file %s",dropFile);
        fprintf(stderr," returned error during write\n");
        exit(1);
      }
      if(fwrite(gain,sizeof(float),NDIM,fpout)<NDIM){
        fprintf(stderr,"lgdtune error: drop file %s",dropFile);
        fprintf(stderr," returned error during write\n");
        exit(1);
      }
      if(fwrite(gc,sizeof(float),NDIM,fpout)<NDIM){
        fprintf(stderr,"lgdtune error: drop file %s",dropFile);
        fprintf(stderr," returned error during write\n");
        exit(1);
      }
      if(fwrite(&nEvents,sizeof(int),1,fpout)<1){
        fprintf(stderr,"lgdtune error: drop file %s",dropFile);
        fprintf(stderr," returned error during write\n");
        exit(1);
      }
      fclose(fpout);
    }

    if(nEvents>1){
      int block;
      gcSolver();    
      for(block=0;block<nChannels;block++){
        if (gc[block]<0.5){
          gain[block]*=0.5;
        }
        else if (gc[block]>2.0){
          gain[block]*=2.0;
        }
        else {
          gain[block]*=gc[block];
        }
      }
      if(outFile){
        FILE *fpout=fopen(outFile,"w");
        if(fpout==NULL){
          fprintf(stderr,"lgdtune error: out file %s",outFile);
          fprintf(stderr," could not be opened for writing\n");
          exit(1);
        }
        for(index=0;index<nChannels;index++){
          row=index/28;
          col=index-(row*28);
          hf1(25,gain[index],1);
          hf2(5,col,row,gain[index]);
	  if(gain[index]!=1){
	    hf2(9,col,row,gain[index]);
	  }
          hf1(6,gc[index],1);
          hf1(7,eigenval[index],1);
          fprintf(fpout," %3d  %2d %2d   %10f   %7f\n",
                  index,col,row,gain[index],gc[index]);
        }
        fclose(fpout);
      }
      fprintf(stderr,"Processed %d events\n",nEvents);
      fprintf(stderr,"Saw data from %d blocks\n",channels);
      fprintf(stderr,"Tuned %d degrees of freedom\n",dof);
    }
  }

/* Create cc output file */

  if(Use_ccOut==1) {
    float dif;
    FILE *fpout=fopen(ccOut,"w");
    for(index=0;index<nChannels;index++){
      cc_out[index] = lgd_cc[index] * gc[index];
      dif = cc_out[index] - lgd_cc[index];
      hf1(50,cc_out[index],1);
      hf1(51,dif,1);
      fprintf(fpout," %f %d\n",cc_out[index],index);
    }
    fprintf(stderr,"Wrote cc info to %s\n",ccOut);
  }
  hrput(0,"lgdtune.hbook","N");
  exit(0);
}
