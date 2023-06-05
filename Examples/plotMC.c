/*
 * dumpMC.c - plot the impacts of the forward particles on the LDG -rtj
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mctypes.h>
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


int main(int argc, char *argv[]);
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
  fprintf(stderr,"  dumpMC <options> file1 file2 ...\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-h\t\tPrint this message\n");
  fprintf(stderr,"\t-M#\t\tProcess at most # events.\n");
}

int main(int argc, char *argv[])

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
	fprintf(stderr,"Reading from %s\n",argv[iarg]);
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
  hrput(0,"plotMC.hbook","N");
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
  tapeHeader_t *hdr=NULL;
  static int currentRun=0;

  if(!event)
    event = malloc(BUFSIZE);
  
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
  mc_vertex_t *mc_vertex=NULL;

  static int hid=0;
  if (hid == 0) {
    hid = 1;
    hbook2(hid, "lgd impacts", 28, -56, 56, 28, -56, 56, 0);
  }
  if (event->eventType == EV_DATA) {
    int i;
    double x, y;
    double z = 118;
    if ((mc_vertex=data_getGroup(event,GROUP_MC_EVENT,1)) == NULL)
      return;
    double momf[3][4];
    for (i=0; i < 2; i++) {
      momf[i][0] = mc_vertex->part[i].momentum.t;
      momf[i][1] = mc_vertex->part[i].momentum.space.x;
      momf[i][2] = mc_vertex->part[i].momentum.space.y;
      momf[i][3] = mc_vertex->part[i].momentum.space.z;
    }
    for (i=0; i < 4; i++) {
      momf[2][i] = momf[0][i] + momf[1][i];
    }
    for (i=0; i < 2; i++) {
      x = momf[i][1] * z / momf[i][3];
      y = momf[i][2] * z / momf[i][3];
      hfill(hid, x, y, 1.);
    }
  }
}
