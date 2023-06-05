/*
  $Log: lgdCalibrate.c,v $
  Revision 1.9  2004/04/02 15:24:51  radphi
  lgdCalibrate.c - removed correctPhotonEnergy() call because makePhotonPatch()
       applies the nonlinearity correction now -rtj-

  Revision 1.8  2001/05/02 22:42:50  radphi
  Changes made by jonesrt@zeus
  minor fixes to allow default for BUFSIZE,
  moved cluster count limit into a define in head of lgdtune.c
  -rtj

  Revision 1.7  2001/04/04 20:22:02  radphi
  Changes made by crsteffe@mantrid00
  updated to use maekHits

  Revision 1.6  2001/03/07 00:40:13  radphi
  Changes made by jonesrt@zeus
  Introduced new lgdtune calibration utility, cleaned up old dat files -rtj

  Revision 1.5  2001/01/11 19:59:30  radphi
  Changes made by armd@urs3
  DSA Jan 11 01    update for using makehits

  Revision 1.4  1999/07/25 03:44:05  radphi
  Changes made by gevmage@jlabs1
  modified to accomodate ccToMap

  Revision 1.3  1999/07/22 04:44:40  radphi
  Changes made by gevmage@jlabs1
  add evt_open and evt_data_read

  Revision 1.2  1998/09/24 16:20:00  radphi
  Changes made by crsteffe@dustbunny
  should work

  Revision 1.1.1.1  1998/06/09 01:02:19  radphi
  Initial revision by radphi@urs2

 * Revision 1.1.1.1  1997/05/24  00:46:54  radphi
 * Initial revision by lfcrob@jlabs2
 *
  */


#include <stdio.h>
#include <stdlib.h>


#include <disData.h>
#include <dataIO.h>
#include <itypes.h>
#include <tapeData.h>
#include <umap.h>
#include <math.h>

#include <lgdGeom.h>
#include <pedestals.h>
#include <eventType.h>
#include <triggerType.h>
#include <cernlib.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <clib.h>
#include <lgdCal.h>
#include <libBGV.h>
#include <makeHits.h>


int main(int argc, char *argv[]);
void Usage(void);
int doAnEvent(int fd);
int processEvent(itape_header_t *event);
int configureIt(itape_header_t *event);
int makeLGDHits(itape_header_t *event);
int doFiducial(lgd_clusters_t *clusters);
void setupHbook(void);
int is_event_clean_of_tdcs(itape_header_t *event, int low_val, int high_val);



#define PI0_MASS 0.135
#define ETA_MASS 0.547

#define TARGET_TO_LGD_DISTANCE 108.0
#define HBOOK 8000000
#define SQR(X) ((X)*(X))
#define E_CRITICAL 0.01455
#define ZMAX_COEFF 4.52
#define ZMAX_MULT .66 /* adjusting the shower depth */
#define SINGLE_BLOCK_CLUSTER_THRESHOLD .150



static float *ped;
static int *thresh;
static float *cc;
static int nChannels,nRows,nCols;
static float baseCC=0.00175;
static int nIterations=1;
static char *ccFile=NULL;

static float massCut=0.0;
static int useMassCut=0;
static float windowCut=0.1;


int ncc[CAL_PRIM_SIZE];
float store[CAL_PRIM_SIZE][CAL_SEC_SIZE];
float x[CAL_PRIM_SIZE];
float u[CAL_PRIM_SIZE];
float newcc[CAL_PRIM_SIZE];

int CCALL_mult_turned_on=0;
float CCALL_mult=1.0;

struct hbk{float ia[HBOOK];} pawc_;


void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  lgdCalibrate <options> file1 file2 ...\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-h\t\tPrint this message\n");
}

int main(int argc, char *argv[])

{
  int iarg;
  char *argptr;
  int fdIn; /* this is the input file descriptor.  evt_open and evt_read_data
		use file descriptors (equivalent to fileno(fdIn) ) instead of 
		the file pointer (fdIn) */
  int maxEvents=0;
  int nEvents=0;
  int index=0;
  int iteration=0;
  int iter2=0;
  FILE *fpCal;
  int icount=0;

  int temprow=0,tempcol=0;

  lgdGeom_t item;

  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      case 'S':  /* Added by Craig Steffen 13MAY98 */
	CCALL_mult=atof(++argptr);
	CCALL_mult_turned_on=1;
	break;
      case 'm':
	massCut = atof(++argptr);
	useMassCut=1;
	break;
      case 'w':
	windowCut = atof(++argptr);
	break;
      case 'c':
	baseCC = atof(++argptr);
	break;
      case 'f':
	ccFile = ++argptr;
	break;
      case 'M':
	maxEvents = atoi(++argptr);
	break;
      case 'h':
	Usage();
	exit(1);
      default:
	fprintf(stderr,"Unknown argument: %s\n",argv[iarg]);
	Usage();
	exit(1);
      }
    }
  ecal0(ncc,store,x,u);
  
  setupHbook();

  for(iteration=0;iteration<nIterations;iteration++){
    for(iarg=1;iarg<argc;iarg++)
      if(*argv[iarg] != '-'){
	if( (fdIn = evt_open(argv[iarg],"r")) == 0){
	  /* 0 is evt_open's error condition */
	  /*	if( (fpIn = fopen(argv[iarg],"r")) == NULL){ */
	  fprintf(stderr,"Failed to open %s for reading\n",argv[iarg]);
	  exit(1);
	}
	else{
	  fprintf(stderr,"Reading from %s\n",argv[iarg]);
	  while( (maxEvents ? maxEvents > nEvents : 1) && !doAnEvent(fdIn)){
	    if( (++nEvents%100) == 0){
	      fprintf(stderr,"%d\r",nEvents);
	      fflush(stderr);
	    }
	  }
	  if(fdIn){
	    evt_close(fdIn);
	    fdIn = 0;
	  }
	}      
      }
    esolv(cc,newcc,ncc,store,x,u);
    fpCal = fopen("cc.dat","w");
    icount=0;
    for(index=0;index<nChannels;index++){
      temprow=index/28;
      tempcol=index-(temprow*28);
      item.channel = index;
      hf2(5,tempcol,temprow,newcc[index]);
      hf1(6,(newcc[index]-cc[index])/cc[index],1);
      hf1(6,newcc[index],1.0);
      fprintf(fpCal," %f %d\n",newcc[index],index);
      hf1(34,newcc[index],1.0);
    }
    fclose(fpCal);
    fprintf(stderr,"Processed %d events\n",nEvents);
    fprintf(stderr,"Calibrated %d modules\n",icount);
  }
  hrput(0,"calibrate.hbook","N");
  return 0;
}

void setupHbook(void)
{
  int index;
  char str[1024];

  hlimit(HBOOK);

  hbook1(1,"Total Energy",200,0.0,8.0,0);
  hbook1(2,"nClusters",10,-0.5,9.5,0);
  hbook2(3,"Hits",28,-0.5,27.5,28,-0.5,27.5,0);  
  hbook2(4,"Hits (Clusterized)",28,-0.5,27.5,28,-0.5,27.5,0);  
  hbook2(5,"New CC",28,-0.5,27.5,28,-0.5,27.5,0);
  hbook1(6,"delta cc",100,0,1,0);

  hbook2(20,"channnel of hit",28,-0.5,27.5,28,-0.5,27.5,0);
  hbook1(21,"radius of clusters (cm)",600,0.0,60.0,0);
  hbook1(22,"radius of photons (cm)",600,0.0,60.0,0);

  hbook1(30,"energy of hit",200,0.0,2.0,0);
  hbook2(31,"E of cluster vs. number of blocks in clusters",13,-0.5,12.5,200,0.0,2.0,0);
  hbook2(32,"position of rejected single block cluster",112,-56.0,56.0,112,-56.0,56.0,0);
  hbook1(33,"E sum of accepted clusters",200,0,5,0);
  hbook1(34,"calibration constants",500,0.0,.005,0);
  hbook1(35,"2 gamma mass for three gamma events,",350,0,1.2,0);

  for(index=1;index<=8;index++){
    /*    sprintf(str,"Total energy, %d clusters",index);
    hbook1(index*100+1,str,200,0,5,0);
    sprintf(str,"Hits, %d clusters",index);
    hbook2(index*100+3,str,28,-0.5,27.5,28,-0.5,27.5,0);
    sprintf(str,"Hits (clusterized), %d clusters",index);
    hbook2(index*100+4,str,28,-0.5,27.5,28,-0.5,27.5,0);
    sprintf(str,"Effective mass, %d clusters",index); */
    hbook1(index*100+20,str,350,0,1.5,0);
    sprintf(str,"Effective mass, fid, %d clusters",index);
    hbook1(index*100+21,str,350,0,1.5,0);
    sprintf(str,"Effective mass, fid, energy, %d clusters",index);
    hbook1(index*100+22,str,350,0,1.5,0);
    sprintf(str,"2 photon effective mass, %d clusters",index);
    hbook1(index*100+30,str,200,0,1.0,0);
  } 
   
  /* Special 2 cluster histos */

  hbook1(2000,"2 gamma chisq",100,0,100,0);

  /* Special 3 cluster histos */

  hbook1(3000,"gamma energy, pi0 gamma",100,0,1,0);
  hbook1(3001,"gamma multipliciy, pi0 gamma",10,0.5,10.5,0);
  hbook1(3002,"gamma pi0 separation, pi0 gamma",80,0,80,0);

  hbook1(4000,"mass, all good 3 cluster",300,0,1.5,0);
  hbook1(4001,"mass, .GE.1 Ok pi0",300,0,1.5,0);
  hbook1(4002,"mass, .EQ.1 good pi0",300,0,1.5,0);
  hbook1(4003,"mass, .EQ.1 good pi0, no good eta",300,0,1.5,0);
  hbook1(4004,"mass, .EQ.1 good pi0, .EQ.1 Ok pi0, no good eta",300,0,1.5,0);
  hbook1(4005,"mass, no good matches",300,0,1.5,0);
  hbook1(4006,"mass, .EQ.1 good eta",300,0,1.5,0);
  hbook1(4007,"mass, .EQ.1 good eta, .EQ.0 good pi0",300,0,1.5,0);
  hbook1(4008,"mass, .EQ.1 good eta, .EQ.0 Ok pi0",300,0,1.5,0);
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
    /* Got an event */
    switch(event->type){
    case TYPE_TAPEHEADER:
      
      /* Get the tape number */
      
      hdr = (void *) event;
      return(doAnEvent(fd));		
    case TYPE_ITAPE:
      if(event->runNo != currentRun){
	/* Load databases */
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
    /* EOF marker on tape */
    fprintf(stderr,"doAnEvent: EOF on tape\n");
    return(doAnEvent(fd));
  case DATAIO_EOT:
    /* Real EOF tape or file */
    return(1);
  case DATAIO_ERROR:
    fprintf(stderr,"doAnEvent: data_read() returned DATAIO_ERROR, Skipping event\n");
    return(1);
  case DATAIO_BADCRC:
    fprintf(stderr,"doAnEvent: data_read() returned DATAIO_BADCRC, Skipping event\n");
    return(1);
  }
  return(1);
}

int configureIt(itape_header_t *event)
{
  int actualRun;
  int index;
  int index2;
  int jpoint;
  float x;

  if(umap_loadDb(event->runNo))
    return(1);
  if(lgdGeomSetup(event->runNo))
    return(1);

  lgdGetGeomParam(&nChannels,&nRows,&nCols);
  

  if(setup_makeHits(event)){
    fprintf(stderr,"Problem with setup_makeHits; exiting");
    return;
  }      
	  
  /* Load pedestals */

  if(ped)
    free(ped);
  ped = malloc(nChannels*sizeof(float));
  if(thresh)
    free(thresh);
  thresh = malloc(nChannels*sizeof(float));

  if(pedestalLoad(PED_LGD,event->runNo,nChannels,ped,thresh,&actualRun))
    return(1);
  printf("pedestalLoad successful for run %d, with %d channels, real run=%d.\n",event->runNo,nChannels,actualRun);
  
  /*  for(index2=0;index2<nChannels;index2++){
    printf("ch=%d,ped=%f,thresh=%d\n",index2,ped[index2],thresh[index2]);
  }*/

  if(!cc){
    cc = malloc(nChannels*sizeof(float));
    if(ccFile){
      FILE *fp;

      fprintf(stderr,"Read cc's from %s\n",ccFile);
      fp = fopen(ccFile,"r");
      for(index=0;index<nChannels;index++){
	cc[index] = baseCC;
      }
      for(index=0;index<nChannels;index++){
	fscanf(fp,"%f %d",&x,&jpoint);
	/* fprintf(stderr," %f %d\n",cc[jpoint],jpoint); */
	cc[jpoint]=x;
      }
      fclose(fp);
      
    }
    else{
      fprintf(stderr,"Set cc's to %f\n",baseCC);
      for(index=0;index<nChannels;index++){
	cc[index] = baseCC;
      }
    }
  }
  /* Added 13MAY98 by Craig Steffen
     change all cc by a constant factor */
  if(CCALL_mult_turned_on){
    for(index=0;index<nChannels;index++){
      cc[index]*=CCALL_mult;
    }
    fprintf(stdout,"configureIt multiplying all cc's by %f.\n",CCALL_mult);
  }

  /* end 13MAY98 change */


  /* Setup the clusterizer */

  /*  lgdSetup(3,0.2,0.1,LGDCLUS_LOG2,0);*/
   lgdSetup(3,0.35,0.15,LGDCLUS_LOG2,0);

  return(0);
}

int processEvent(itape_header_t *event)
{

  lgd_hits_t *lgdHits=NULL;
  float Etot=0.0,Etot_clus=0.0;
  int index,index2;
  lgd_clusters_t *clusters=NULL,*new_clusters=NULL;
  lgd_hits_t *clusterHits=NULL,*new_clusterHits=NULL;
  void *tmp=NULL;
  int size;
  vector4_t pGamma[8];
  vector3_t vertex={0.0,0.0,0.0};
  float mass;
  float mass2;
  int calMassCut=0;

  int N_good_pi0=0;
  int N_Ok_pi0=0;
  int N_good_eta=0;
  int iter1;

  int single_block_cluster_flag=0;

  int temprow,tempcol;

  float sigma1,sigma2;
  vector4_t pFit[8];
  lgdGeom_t item;
  float radiusMin;
  float radius;
  float mag;
  int k;
  int threeCombo[3][2]={{1,2},{0,2},{0,1}};
  int cluster_iter,hit_iter;

  if( (event->eventType == EV_DATA) && (event->trigger == TRIG_DATA)){
    if(unpackEvent(event,BUFSIZE))
       return(1);
    /* this eliminates TDC hits */
    if(!is_event_clean_of_tdcs(event,420,480))
      return(1);
    /*    if(makeLGDHits(event))
	  return(1); */
    if(makeHits(event)){
      fprintf(stderr,"makeHits returned an error!\n");
      return(1);
    }
    lgdHits = data_getGroup(event,GROUP_LGD_HITS,0);
    /*    printf("there were %d hits.\n",lgdHits->nhits); */
    for(index=0;index<lgdHits->nhits;index++){
      Etot += lgdHits->hit[index].energy;
      temprow=lgdHits->hit[index].channel/28;
      tempcol=lgdHits->hit[index].channel-(temprow*28);
      hf2(20,(float)temprow,(float)tempcol,1.0);
    }
    hf1(1,Etot,1);

    /* Clusterize */

    lgdClusterizeIU(lgdHits,&clusterHits,&clusters);

    /* Craig Steffen 19MAY98
       removing one block clusters with < .25 of energy */
    
    new_clusterHits=(lgd_hits_t *)malloc(sizeof(lgd_hits_t)+((clusterHits->nhits)*sizeof(lgd_hit_t)));
    new_clusters=(lgd_clusters_t *)malloc(sizeof(lgd_clusters_t)+((clusters->nClusters)*sizeof(lgd_cluster_t)));
    new_clusterHits->nhits=0;
    new_clusters->nClusters=0;
    for(cluster_iter=0;cluster_iter<clusters->nClusters;cluster_iter++){
      if((clusters->cluster[cluster_iter].nBlocks==1)&&(clusters->cluster[cluster_iter].energy<SINGLE_BLOCK_CLUSTER_THRESHOLD)){
	/* rejected cluster stuff here */
	hf2(32,clusters->cluster[cluster_iter].space.x,clusters->cluster[cluster_iter].space.y,1.0);
      }
      else{
	new_clusters->cluster[new_clusters->nClusters].nBlocks=clusters->cluster[cluster_iter].nBlocks;
	new_clusters->cluster[new_clusters->nClusters].flags=clusters->cluster[cluster_iter].flags;
	new_clusters->cluster[new_clusters->nClusters].space=clusters->cluster[cluster_iter].space;
	new_clusters->cluster[new_clusters->nClusters].width=clusters->cluster[cluster_iter].width;
	new_clusters->cluster[new_clusters->nClusters].energy=clusters->cluster[cluster_iter].energy;
	new_clusters->cluster[new_clusters->nClusters].firstBlock=new_clusterHits->nhits;
	for(hit_iter=(clusters->cluster[cluster_iter].firstBlock);hit_iter<((clusters->cluster[cluster_iter].firstBlock)+(clusters->cluster[cluster_iter].nBlocks));hit_iter++){
	  new_clusterHits->hit[new_clusterHits->nhits].channel=clusterHits->hit[hit_iter].channel;
	  new_clusterHits->hit[new_clusterHits->nhits].energy=clusterHits->hit[hit_iter].energy;
	  (new_clusterHits->nhits)++;
	}
	(new_clusters->nClusters)++;
      }
    }
    free(clusterHits);
    free(clusters);
    clusterHits=NULL;
    clusters=NULL;
    clusterHits=new_clusterHits;
    clusters=new_clusters;

    data_removeGroup(event,GROUP_LGD_CLUSTERS,0);
    data_removeGroup(event,GROUP_LGD_CLUSTER_HITS,0);
    
    size = (int)sizeof(lgd_clusters_t) + (clusters->nClusters-1)*(int)sizeof(lgd_cluster_t);
    tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_CLUSTERS,0,size);
    if(!tmp){
      fprintf(stderr,"processEvent: Failed to add clusters\n");
      free(clusters);
      free(clusterHits);
      return(1);
    }
    memcpy(tmp,clusters,size);
    
    tmp=NULL;
    size = (int)sizeof(lgd_hits_t) + (clusterHits->nhits)*(int)sizeof(lgd_hit_t);
    tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_CLUSTER_HITS,0,size);
    if(!tmp){
      fprintf(stderr,"processEvent: Failed to add cluster hits\n");
      free(clusters);
      free(clusterHits);
      return(1);        
    }
    memcpy(tmp,clusterHits,size);
    /*    if(clusters->nClusters)
	  fprintf(stderr,"%d Clusters\n",clusters->nClusters);
      */
    fflush(stdout);

    hf1(2,clusters->nClusters,1);
    /*    hf1(100*clusters->nClusters+1,Etot,1); */
    for(index=0;index<lgdHits->nhits;index++){
      item.channel = lgdHits->hit[index].channel;
      if(!lgdGetGeom(0,&item,LGDGEOM_BYCHANNEL)){
	hf2(3,item.col,item.row,1);
	/*	hf2(100*clusters->nClusters+3,item.col,item.row,1); */
      }
    }
    Etot=0.;
    for(index=0;index<clusters->nClusters;index++){
      /*
            pGamma[index]=makePhotonGen(MAKE_PHOTON_DEPTH,&clusters->cluster[index],NULL,vertex); 
      */
      /* 
       * Scotts Kludge !!!!
       */
	    /*                       clusters->cluster[index].space.z=128.;
      mag=sqrt(clusters->cluster[index].space.x*clusters->cluster[index].space.x
	       +clusters->cluster[index].space.y*clusters->cluster[index].space.y
	       +clusters->cluster[index].space.z*clusters->cluster[index].space.z);
      pGamma[index].t=1.23888*pGamma[index].t;
      Etot+=pGamma[index].t;
      pGamma[index].space.x=clusters->cluster[index].space.x*pGamma[index].t/mag;
      pGamma[index].space.y=clusters->cluster[index].space.y*pGamma[index].t/mag;
      pGamma[index].space.z=clusters->cluster[index].space.z*pGamma[index].t/mag; */
      /*
       * end of Kluge
       */
      
      /* Correct algorithm from E852 code.  Needs to be placed into 
	 libraries and stuff.  Craig Steffen 11MAY98 */

      /*      printf("%f\t%f\t%f\t%f\n",clusters->cluster[index].energy,clusters->cluster[index].space.x,clusters->cluster[index].space.y,clusters->cluster[index].space.z);
      printf("%f\t%f\t%f\t%f\n",pGamma[index].t,pGamma[index].space.x,pGamma[index].space.y,pGamma[index].space.z);
      */

          /* do depth correction; use library routine  */
      pGamma[index]=makePhotonPatch(&(clusters->cluster[index]),vertex);

          /* do non-linearity correction; use library routine  */

 /* non-linearity is now treated together with depth correction
  * in makePhotonPatch() above -rtj-
      pGamma[index] = correctPhotonEnergy(pGamma[index]); 
  */

      Etot+=pGamma[index].t;

      
      
      /*      printf("%f\t%f\t%f\t%f\n",pGamma[index].t,pGamma[index].space.x,pGamma[index].space.y,pGamma[index].space.z); */
    }
    hf1(33,Etot,1.0);
    fflush(stdout);
    /*        printf("\n"); */
    for(index=0;index<clusters->nClusters;index++){
      for(index2=index+1;index2<clusters->nClusters;index2++){
	mass = nPartMass(2,&pGamma[index],&pGamma[index2]);
	hf1(100*clusters->nClusters+30,mass,1);
      }
      hf2(31,clusters->cluster[index].nBlocks,clusters->cluster[index].energy,1.0);
    }
    
    for(index=0;index<clusterHits->nhits;index++){
      item.channel = clusterHits->hit[index].channel;
      if(!lgdGetGeom(0,&item,LGDGEOM_BYCHANNEL)){
	hf2(4,item.col,item.row,1);
	/*	hf2(100*clusters->nClusters+4,item.col,item.row,1); */
      }
    }
    /*    printf("+");
    fflush(stdout);
    */
    switch(clusters->nClusters){
    case 2:
      /*    printf("2"); */
    fflush(stdout);
      mass = nPartMass(2,&pGamma[0],&pGamma[1]);
      hf1(220,mass,1);
      hf1(230,mass,1);
      if(!doFiducial(clusters)){
	/*printf("f");*/
	fflush(stdout);
	hf1(221,mass,1);
	if(Etot> 1.0){
	  hf1(222,mass,1);
	  calMassCut=0;
	  if(useMassCut){
	    if(mass < massCut){
	      calMassCut = 1;
	    }
	  }
	  else{
	    if(fabs(mass-PI0_MASS) < windowCut){
	      calMassCut=1;
	    }
	  }
	  if(calMassCut){
	    float chisq;
	    sigma1 = 0.02*pGamma[0].t + 0.08*sqrt(pGamma[0].t);
	    sigma2 = 0.02*pGamma[1].t + 0.08*sqrt(pGamma[1].t);
	    chisq = fitTwoGamma(&pGamma[0],&pGamma[1],PI0_MASS,
				sigma1,sigma2,&pFit[0],&pFit[1]);
	    hf1(2000,chisq,1);
	    if(chisq > 50)
	      return(1);
	    if( (pFit[0].t < 0) || (fabs(pFit[0].t) < 0.00001) || (is_NaN(pFit[0].t))){
	      fprintf(stderr,"Bad fitted energy 0 (%f)!\n",pFit[0].t);
	      return(1);
	    }
	    if( (pFit[1].t < 0) || (fabs(pFit[1].t) < 0.00001) || (is_NaN(pFit[1].t))){
	      fprintf(stderr,"Bad fitted energy 1 (%f)!\n",pFit[1].t);
	      return(1);
	    }
	    ecal(pFit[0].t,&clusters->cluster[0],
		 clusterHits,cc,ncc,store,x,u);
	    ecal(pFit[1].t,&clusters->cluster[1],
		 clusterHits,cc,ncc,store,x,u);
	    
	  }
	}
      }
      /*      printf(")"); */
      fflush(stdout);
      break;
    case 3:
      /*    printf("3"); */
    fflush(stdout);
      mass = nPartMass(3,&pGamma[0],&pGamma[1],&pGamma[2]);
      /*      printf("\nM_t=%f ",mass); */
      hf1(320,mass,1);
      if(!doFiducial(clusters)){
	/*	printf("\n");
	for(iter1=0;iter1<3;iter1++){
	  printf("E=%f  p_x=%f  p_y=%f  p_z=%f\n",pGamma[iter1].t,pGamma[iter1].space.x,pGamma[iter1].space.y,pGamma[iter1].space.z);
	} */
	/* omega search */
	hf1(321,mass,1);
        if(Etot>2.0){
	  hf1(322,mass,1);
	}
	for(k=0;k<3;k++){
	  mass2= nPartMass(2,&pGamma[threeCombo[k][0]],&pGamma[threeCombo[k][1]]);
	  /*	  printf("M_%d=%f ",k,mass2); */
	  /* WARNING! HARD CODED NUMBER!!! Craig Steffen */
	  /* changed from .04 to .02 21MAY98 */
	  /* changed to .04 again 17JUL98 */
	  if((fabs(mass2-PI0_MASS) < 0.04) && (clusters->cluster[k].nBlocks != 1)){
	    hf1(3000,pGamma[k].t,1);
	    hf1(3001,clusters->cluster[k].nBlocks,1);
	    
	    radiusMin=99999.0;
	    for(index=0;index<3;index++){
	      if(index!=k){
		radius = SQR(clusters->cluster[index].space.x - clusters->cluster[k].space.x) + SQR(clusters->cluster[index].space.y - clusters->cluster[k].space.y);
		
		if(radius < radiusMin)
		  radiusMin= radius;	      
	      }	      
	    }	
	    hf1(3002,sqrt(radiusMin),1);
	    /*	    if(Etot>3.0){
	      hf1(3003,mass,1);
	    } */
	  } /* if((fabs(mass2-PI0_MASS) < 0.01) && (clusters->cluster[k].nBlocks != 1)) */
	  if((fabs(mass2-ETA_MASS) < 0.1) && (clusters->cluster[k].nBlocks != 1)){
	    /*	    if(Etot>3.0){
	      hf1(3013,mass,1);
	    } */
	  } /*  eta mass check */
	  if((Etot>3.0) && (clusters->cluster[k].nBlocks != 1)){
	    if(fabs(mass2-PI0_MASS) < 0.08){
	      N_Ok_pi0++;
	      if(fabs(mass2-PI0_MASS) < 0.04){
		N_good_pi0++;
	      }
	    }
	    if(fabs(mass2-ETA_MASS) < 0.10){
	      N_good_eta++;
	    }
	  }
	} /* loops over three combinations */
	
	/* Craig's anti-combinatoric code 23AUG98 */
	hf1(4000,mass,1.0);
	if(N_Ok_pi0>=1){
	  hf1(4001,mass,1.0);
	}
	if(N_good_pi0==1){
	  hf1(4002,mass,1.0);
	  if(N_good_eta==0){
	    hf1(4003,mass,1.0);
	    if(N_Ok_pi0==1){
	      hf1(4004,mass,1.0);
	    }
	  }
	}
	if(N_good_eta==1){
	  hf1(4006,mass,1.0);
	  if(N_good_pi0==0){
	    hf1(4007,mass,1.0);
	    if(N_Ok_pi0==0){
	      hf1(4008,mass,1.0);
	    }
	  }
	}
	if(!(N_good_eta+N_Ok_pi0)){
	  hf1(4005,mass,1.0);
	}
      } /* if passes fiducial cut */
      break;
    case 4:
      /*    printf("4"); */
      fflush(stdout);
      mass = nPartMass(4,&pGamma[0],&pGamma[1],&pGamma[2],&pGamma[3]);
      hf1(420,mass,1);
      if(!doFiducial(clusters)){
	hf1(421,mass,1);
      }
    default:
      break;
    }
    /*    printf("d"); */
    fflush(stdout);
    
    
    free(clusters);
    free(clusterHits);
  }
  /*  else{
    if(event->eventType == EV_DATA)
      fprintf(stderr,"Non-physics event: 0x%x\n",event->trigger);
  }*/
  
  /*  printf("@"); */
  fflush(stdout);
}

int doFiducial(lgd_clusters_t *clusters)
{
  int index;
  float radius;

  for(index=0;index<clusters->nClusters;index++){
    radius = SQR(clusters->cluster[index].space.x) + SQR(clusters->cluster[index].space.y);
    if(radius < 100)
      return(1);
  }
  return(0);
  
}

int makeLGDHits(itape_header_t *event)
{
  adc_values_t *lgd=NULL;
  static lgd_hits_t *lgdHits=NULL;
  static lgd_hits_t *tmp=NULL;
  int index;
  int size;
  lgdGeom_t item;
  
  if(!lgdHits){
    lgdHits = malloc(sizeof(adc_hits_t)*(nChannels-1)*sizeof(adc_hit_t));
  }

  lgdHits->nhits = 0;

  if( (lgd = data_getGroup(event,GROUP_LGD_ADCS,0)) == NULL){
    fprintf(stderr,"No LGD ADC's\n");
    return(1);
  }
  for(index=0;index<lgd->nadc;index++){
    /*    printf("index=%d,value=%d, thresh=%d, pedestal=%f\n",index,lgd->adc[index].value,thresh[lgd->adc[index].channel],ped[lgd->adc[index].channel]);*/
    if(lgd->adc[index].value > thresh[lgd->adc[index].channel]){ 
    /*    if((((float)lgd->adc[index].value) - ped[lgd->adc[index].channel])>8.0 ){ */
      lgdHits->hit[lgdHits->nhits].channel = lgd->adc[index].channel;
      lgdHits->hit[lgdHits->nhits].energy = (lgd->adc[index].value - ped[lgd->adc[index].channel])*cc[lgd->adc[index].channel];      
      hf1(30,lgdHits->hit[lgdHits->nhits].energy,1.0);
      lgdHits->nhits++;
    }
  }
  /*  printf("nhits=%d\n",lgdHits->nhits);*/
  if(lgdHits->nhits)
    size = sizeof(lgd_hits_t) + (lgdHits->nhits-1)*sizeof(lgd_hit_t);
  else
    size = sizeof(lgd_hits_t);
  tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_HITS,0,size);
  memcpy(tmp,lgdHits,size);
  return(0);
}

int is_event_clean_of_tdcs(itape_header_t *event, int low_val, int high_val){
  /* must be given unpacked data */
  tdc_values_t *cpv_tdcs=NULL;
  int i;
  cpv_tdcs=data_getGroup(event,GROUP_CPV_TDCS,0);
  if(cpv_tdcs==NULL)
    return 1; /* if on TDC group, assume is clear */
  for(i=0;i<cpv_tdcs->ntdc;i++){
    if(cpv_tdcs->tdc[i].le>low_val&&cpv_tdcs->tdc[i].le<high_val)
      return 0;
  }
  return 1;
}

