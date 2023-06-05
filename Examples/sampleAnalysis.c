/*

  $Log: sampleAnalysis.c,v $
  Revision 1.7  2014/06/06 21:51:05  radphi
  * merge_clusters.c [rtj]
     - added new function to allow merging of low-mass clusters in the LGD.
  * ntuple.c [rtj]
     - implement low-mass cluster merging in the LGD if the number of photons
       is larger than 3.

  Revision 1.6  2001/03/07 00:51:43  radphi
  Changes made by jonesrt@zeus
  Additions to mctuple and ntuple to include new mesons group -rtj

  Revision 1.5  2000/01/18 03:52:41  radphi
  Changes made by jonesrt@golem
  ported Examples codes: sampleAnalysis, dumpEvent and
  paw helper routines scalers.sl, mctuple.sl, ntuple.sl to linux-Redhat (Intel)
  -rtj

  Revision 1.4  1999/06/23 23:04:03  radphi
  Changes made by jonesrt@jlabs1
  converted sizeof to use macros - RTJ,TBB

 * Revision 1.3  1998/04/17  23:55:58  radphi
 * Changes made by jonesrt@jlabs3
 * just bringing scalers.c up to date - R.T.Jones
 *
 * Revision 1.2  1997/07/28  12:03:09  radphi
 * Changes made by lfcrob@jlabs2
 * Fixed uninitialized pointers
 *
 * Revision 1.1  1997/06/02  21:07:10  radphi
 * Initial revision by lfcrob@dustbunny
 * Kind of messy shell for an analysis program
 *

  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <disData.h>
#include <dataIO.h>
#include <itypes.h>
#include <tapeData.h>
#include <umap.h>
#include <math.h>
#include <unpackData.h>

#include <lgdGeom.h>
#include <cpvGeom.h>
#include <pedestals.h>
#include <calibration.h>
#include <eventType.h>
#include <triggerType.h>
#include <cernlib.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <clib.h>


void main(int argc, char *argv[]);
void Usage(void);
int doAnEvent(FILE *fp);
void processEvent(itape_header_t *event);
int configureIt(itape_header_t *event);
int makeLGDHits(itape_header_t *event);
int doFiducial(lgd_clusters_t *clusters);
void setupHbook(void);

#define BUFSIZE 100000
#define HBOOK 8000000
#define SQR(X) ((X)*(X))

static float *pedLGD=NULL;
static float *pedUPV=NULL;
static int *threshLGD=NULL;
static int *threshUPV=NULL;
static float *cc=NULL;
static int nChannelsUPV,nE_UPV,nG_UPV;
static int nChannelsLGD,nRowsLGD,nColsLGD;
static int nChannelsCPV,nHCPV,nVCPV;


/* For PAW */

struct hbk{float ia[HBOOK];} pawc_;


void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  sampleAnalysis <options> file1 file2 ...\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-h\t\tPrint this message\n");
  fprintf(stderr,"\t-M#\t\tProcess at most # events.\n");
}

void main(int argc, char *argv[])

{
  int iarg;
  char *argptr;
  FILE *fpIn;
  int maxEvents=0;
  int nEvents=0;
  int index=0;
  int iteration=0;
  FILE *fpCal;

  lgdGeom_t item;

  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
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
  
  setupHbook();

  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] != '-'){
      if( (fpIn = fopen(argv[iarg],"r")) == NULL){
	fprintf(stderr,"Failed to open %s for reading\n",argv[iarg]);
	exit(1);
      }
      else{
	fprintf(stderr,"Reading from %s (fileno %d)\n",argv[iarg],fileno(fpIn));
	while( (maxEvents ? maxEvents > nEvents : 1) && !doAnEvent(fpIn)){
	  if( (++nEvents%1000) == 0){
	    fprintf(stderr,"%d\r",nEvents);
	    fflush(stderr);
	  }
	}
	if(fpIn){
	  fclose(fpIn);
	  fpIn = NULL;
	}
      }      
    }  
  fprintf(stderr,"Processed %d events\n",nEvents);
  hrput(0,"sample.hst","N");
}

void setupHbook(void)
{
  hlimit(HBOOK);

  /* Put hbook() calls here */

}
int doAnEvent(FILE *fp)
{
  static itape_header_t *event=NULL;
  int ret;
  int fd;
  tapeHeader_t *hdr=NULL;
  static int currentRun=0;

  if(!event)
    event = malloc(BUFSIZE);
  
  fd = fileno(fp);
  if (fd == -1) {
    fprintf(stderr,"Input file pointer fp is invalid!\n");
    return(1);
  }
  ret=data_read(fileno(fp),event,BUFSIZE);
  switch(ret){      
  case DATAIO_OK:
    /* Got an event */
    switch(event->type){
    case TYPE_TAPEHEADER:
      
      /* Get the tape number */
      
      hdr = (void *) event;
      fprintf(stderr,"Prcoessing tape %d\n",hdr->tapeSeqNumber);
      return(doAnEvent(fp));		
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
      return(doAnEvent(fp));
    }
  case DATAIO_EOF:
    /* EOF marker on tape */
    fprintf(stderr,"doAnEvent: EOF on tape\n");
    return(doAnEvent(fp));
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

  if(umap_loadDb(event->runNo))
    return(1);
  if(lgdGeomSetup(event->runNo))
    return(1);
  if(cpvGeomSetup(event->runNo))
    return(1);

  lgdGetGeomParam(&nChannelsLGD,&nRowsLGD,&nColsLGD);
  cpvGetGeomParam(&nChannelsCPV,&nHCPV,&nVCPV);
  
  /* Load pedestals */

  if(pedLGD)
    free(pedLGD);
  pedLGD = malloc(nChannelsLGD*sizeof(float));
  if(threshLGD)
    free(threshLGD);
  threshLGD = malloc(nChannelsLGD*sizeof(float));
  if(pedestalLoad(PED_LGD,event->runNo,nChannelsLGD,pedLGD,threshLGD,&actualRun))
    return(1);


  if(cc)
    free(cc);
  cc = malloc(nChannelsLGD*sizeof(float));
  if(ccLoad(CC_LGD,event->runNo,nChannelsLGD,cc,&actualRun))
    return(1);

  /* Setup the clusterizer */

  lgdSetup(3,0.2,0.05,LGDCLUS_LOG2,0);

  return(0);
}

void processEvent(itape_header_t *event)
{

  lgd_hits_t *lgdHits=NULL;
  adc_values_t *lgd_adcs;
  tdc_values_t *tagger;
  float Etot=0;
  int index,index2;
  lgd_clusters_t *clusters=NULL;
  lgd_hits_t *clusterHits=NULL;
  void *tmp=NULL;
  int size;
  vector4_t pGamma[8];
  vector3_t vertex={0,0,0};
  float mass;
  float mass2;
  int calMassCut=0;

  float sigma1,sigma2;
  vector4_t pFit[8];
  lgdGeom_t item;
  float radiusMin;
  float radius;

  int k;
  int threeCombo[3][2]={{2,3},{1,3},{1,2}};

  if( (event->eventType == EV_DATA) && (event->trigger == TRIG_DATA)){
    if(unpackEvent(event,BUFSIZE))
      return;

    if((lgd_adcs=data_getGroup(event,GROUP_LGD_ADCS,0)) == NULL)
      return;

/*
  This little section will print out an ASCII map of the LGD
  as seen standing behind the wall and looking upstream, with
  the counters that were mapped once to row,col addresses
  showing their correct channel count, the missing blocks
  showing up with little dots, and the ones mapped multiple
  times shown with asterixes.             R.T.Jones  3-Apr-98
*/
    if (0) {
      int map[800];
      int i, j, chan;
      for (chan=0; chan < 800; chan++)
        map[chan] = 0;
      for (i=0; i < lgd_adcs->nadc; i++)
        map[lgd_adcs->adc[i].channel]++;
      for (i=27; i > -1; i--) {
        for (j=0; j < 28; j++) {
          chan = i*28+j;
          switch (map[chan]) {
          case 0:
            printf(" .. ");
            break;
          case 1:
            printf("%3d ",chan);
            break;
          default:
            printf(" ** ");
          }
        }
        printf("\n\n");
      }
    }

    if(makeLGDHits(event))
      return;

    /* Get tagger TRC's */

    if((tagger=data_getGroup(event,GROUP_TAGGER_TDCS,0)) == NULL)
      return;

    lgdHits = data_getGroup(event,GROUP_LGD_HITS,0);

    /* Clusterize */

    lgdClusterizeIU(lgdHits,&clusterHits,&clusters);
    data_removeGroup(event,GROUP_LGD_CLUSTERS,0);
    data_removeGroup(event,GROUP_LGD_CLUSTER_HITS,0);
    
    size = (int)sizeof_lgd_clusters_t(clusters->nClusters);
    tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_CLUSTERS,0,size);
    if(!tmp){
      fprintf(stderr,"processEvent: Failed to add clusters\n");
      free(clusters);
      free(clusterHits);
      return;
    }
    memcpy(tmp,clusters,size);
    
    tmp=NULL;
    size = (int)sizeof_lgd_hits_t(clusterHits->nhits);
    tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_CLUSTER_HITS,0,size);
    if(!tmp){
      fprintf(stderr,"processEvent: Failed to add cluster hits\n");
      free(clusters);
      free(clusterHits);
      return;        
    }
    memcpy(tmp,clusterHits,size);


    /* Make photons */

    for(index=0;index<clusters->nClusters;index++){
      pGamma[index]=makePhotonGen(MAKE_PHOTON_DEPTH,&clusters->cluster[index],NULL,vertex);
    }

    /* Calculate 2 photon mass */

    if(clusters->nClusters==2){
      mass = nPartMass(2,&pGamma[0],&pGamma[1]);
    }    
    free(clusters);
    free(clusterHits);
  }
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
    lgdHits = malloc(sizeof_adc_hits_t(nChannelsLGD));
  }

  lgdHits->nhits = 0;

  if( (lgd = data_getGroup(event,GROUP_LGD_ADCS,0)) == NULL){
    fprintf(stderr,"No LGD ADC's\n");
    return(1);
  }
  for(index=0;index<lgd->nadc;index++){
    if(lgd->adc[index].value > threshLGD[lgd->adc[index].channel]){
      lgdHits->hit[lgdHits->nhits].channel = lgd->adc[index].channel;
      lgdHits->hit[lgdHits->nhits].energy = (lgd->adc[index].value - pedLGD[lgd->adc[index].channel])*cc[lgd->adc[index].channel];      
      lgdHits->nhits++;
    }
  }
  if(lgdHits->nhits)
    size = sizeof_lgd_hits_t(lgdHits->nhits);
  else
    size = sizeof_lgd_hits_t(1);
  tmp = data_addGroup(event,BUFSIZE,GROUP_LGD_HITS,0,size);
  memcpy(tmp,lgdHits,size);
  return(0);
}





