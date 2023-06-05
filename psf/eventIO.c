#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <dataIO.h>
#include <disData.h>
#include <makeHits.h>
#include <triggerType.h>
#include <detector_params.h>

#include <ntypes.h>
#include <mctypes.h>
#include "../Examples/mctuple.h"

#include <eventIO.h>
#include <fileIO.h>
#include <cernlib.h>

/* private data members */

#define FILE_LIST "file.in"

#define SHOWER_IMAGE_DIM 9
#define SHOWER_IMAGE_OFFSET 4
#define SHOWER_CENTER_NEIGHBOUR 2
#define LGD_BLOCK_SIZE 4
#define LGD_MATRIX_ROWS 28
#define LGD_MATRIX_COLS 28
#define LGD_PLANE_DISTANCE 120
#define THE_PHOTON_PHI 0.2
#define MC_RUN 9999
#define SQR(X) ((X)*(X))

//#define ISOLATED_SHOWERS 
//#define GET_MC_FROM_FILE_NAME
//#define VERBOSE

int fd=0;
int runNumber=0;
int eventNumber=0;
int iPhoton=-1;
int iPhotons=0;
int ndim = LGD_MATRIX_ROWS*LGD_MATRIX_COLS;
int isMC=0;
int isLgdImg=0;
int isShowerImg=0;
int isFirst=1;
int define=0;

FILE *fList=NULL;
char *inFileName=NULL;
char *mcFile=NULL;

lgd_hits_t *lgdHits=NULL;
itape_header_t *event=NULL;
photons_t *photons=NULL;
mc_event_t *mc_event=NULL;
mc_vertex_t *mc_vertex=NULL;
nt_mc_event_t *p_mc=NULL;

struct {
   int centerChan;		/* central channel of the image */
   int centerRow;		
   int centerCol;     
   float Emax;                  /* energy of maximum block  */
   float Etot;			/* energy sum in the image */
   float yield[SHOWER_IMAGE_DIM][SHOWER_IMAGE_DIM];
} lgd_shower_image;

struct { 
  int maxChan; 			/* max energy block in the LGD*/
  int maxRow;
  int maxCol;
  float Emax;
  float Etot;			/* energy sum in the LGD */
  float yield[LGD_MATRIX_ROWS][LGD_MATRIX_COLS]; 
} lgd_image;


/* public function members defined in eventIO.h */

/*  public function members (c interface) 
int eventGetNext(void);
* get LGD characteristic variables 
float eventGetLgdEtot(void);
float eventGetLgdEmax(int *ch);
float eventGetLgdE(int ch);
float eventGetLgdErc(int row, int col);
* get shower characteristic variables 
float eventGetEtot(void);
float eventGetEmax(int *ch);
float eventGetE(int ch);
float eventGetErc(int row, int col);
* get X,Y and R,C from LGD block (channel) 
void eventGetLgdXY(int row, int col, float *x, float *y);
void eventGetLgdRC(int ch, int *row, int *col);

* public function members (fortran interface) 
int eventtest_(int *count);
void eventclose_(void);
int eventnext_(void);
int eventmc_(void);
void eventimage_(void);
float eventphoton_(float *Px, float *Py, float *Pz);
float mcphoton_(float *Px, float *Py, float *Pz);
int eventrun_();
*/

/* private function members */
static int getNextFileIn(void);
static lgd_hits_t *eventSetup(itape_header_t *event);
static photons_t *eventPhotons(itape_header_t *event);
//static mc_event_t *eventMCdata(itape_header_t *event);
static int eventIsMC(itape_header_t *event);
static int eventSetP_MC(itape_header_t *event);
static int eventGetPhotonAngle4Name(float *theta, float *phi);
static float eventGetPhotonEmc4Name(float *theta, float *phi);
static int eventFindPhoton(float theta, float phi);
static float eventGetPhotonE(int *list, float *Px, float *Py, float *Pz);
static float eventGetPhotonEmc(int *list, float *Px, float *Py, float *Pz);
static int makeLgdImage(void);
static int makeShowerImageHot(void);
static int makeShowerImagePhot(int iphot);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * 
* Private method to obtain the name of the file
* from the list defined by FILE_LIST (file.in by default)
*/
static int getNextFileIn(void) {
  int ret;
  if ((fList == NULL) && (fList=fopen(FILE_LIST,"r")) == NULL) {
    fprintf(stderr,"ERROR - cannot open file list \n");
    return -1;
  }
  if (inFileName == NULL) {
    inFileName = malloc(1000);
  }
  if ( ret=fscanf(fList,"%s",inFileName) == EOF ) {
    fprintf(stderr,"End of input file list\n");
    free(inFileName);
    inFileName = NULL;
    fclose(fList);
    fList=NULL;
    return 0;
  }
  return 1;
}

/* * * * * * * * * * * * * * * * * * * * * * * 
* Private method to setup an event and to 
* return a pointer to lgdHits
*
* -revision: 8.6.02 - in addition get clusters and photons groups
*		      and put photons in the eventIO global space. 
*/
static lgd_hits_t *eventSetup(itape_header_t *event) {

  if ( (event->eventType) != EV_DATA ) {
     fprintf(stderr,"setupEvent: exit, not an EV_DATA! \n");
     return NULL;
  }

  if ( event->runNo != runNumber ) {
    if(setup_makeHits(event)){
      fprintf(stderr,"ERROR - problem with setup_makeHits\n");
      return NULL;
    }
    runNumber = event->runNo;
  }

/* unpack event just for real data */
  if ( event->trigger == TRIG_DATA) {
    if(unpackEvent(event,BUFSIZE)){
      fprintf(stderr,"ERROR -  unpacking error!\n");
      return NULL;
    }
    
    if(make_lgd_hits(event,NULL,640)){
#ifdef VERBOSE
      fprintf(stderr,"ERROR -  make_lgd_hits error!\n");
#endif
      return NULL;
    }
  }
  
  return data_getGroup(event,GROUP_LGD_HITS,0);
}


/* * * * * * * * * * * * * * * * 
* Private method to set a pointer to the Photons Group
*/
static photons_t *eventPhotons(itape_header_t *event) {

  if(make_lgd_clusters(event,NULL,640)) {
#ifdef VERBOSE
    fprintf(stderr,"ERROR -  make_lgd_clusters error!\n");
#endif
    return NULL;
  }
  if(make_photons(event,NULL,640)) {
#ifdef VERBOSE
    fprintf(stderr,"ERROR -  make_photons error!\n");
#endif
    return NULL;
  }
  
  return data_getGroup(event,GROUP_PHOTONS,0);
}



/* * * * * * * * * * * * * * * * 
* Private method to determine if this is a Monte Carlo event and set
* pointer to MC data group.
*  - in case of simulation that contains gradphi kine data embeded 
* in the itape file name set the pointer to the file name and return 3.
*  - in the case of old single shower generation that does not have 
* MC group return 2. 
*  - ow return 1.
*/
static int eventIsMC(itape_header_t *event) {
    
  if (event->runNo != MC_RUN) {
    return 0;
  }
#ifdef GET_MC_FROM_FILE_NAME
  if (mcFile == NULL) {
    mcFile = malloc(1000);
  }
  if ( strcmp(mcFile,inFileName) ) {
    strcpy(mcFile,inFileName);
  }
  return 3;
#endif
  if ((mc_event=data_getGroup(event,GROUP_MC_EVENT,0)) == NULL) {
    return 2;  
  }
  return 1;
}
/* * * * * * * * * * * * * * * * 
* Private method to obtain MC group from event
*/  
static mc_event_t *eventMCdata(itape_header_t *event) {
  
  mc_event_t *pMC=NULL;

  if(isMC == 0) { 
    fprintf(stderr,"Warning - no MC data");	  
    return NULL;
  }
  
  if ((pMC=data_getGroup(event,GROUP_MC_EVENT,0)) == NULL) {
    fprintf(stderr,"ERROR - failed in geting MC data");
    return NULL;
  }
  return pMC;
}

/* * * * * * * * * * * * * * * * 
* Private method to set a pointer to the MC data structure
* like the one used in ntuple (mctuple).
* New simulation production should contain MC group!
* Returns 1 on success, otherwise zero!
*/ 
static int eventSetP_MC(itape_header_t *event){

 int i,j;

 if (mc_event == NULL) {return -1;}
 
 if (p_mc)    free(p_mc), p_mc = NULL;
 
 if ((p_mc=malloc(sizeof(nt_mc_event_t))) == NULL) {
   fprintf(stderr,"ERROR - could not alocate p_mc\n");
   return -1;
 }  
   
 p_mc->nI=mc_event->npart;
 p_mc->type=mc_event->type;
 p_mc->wgt=mc_event->weight;
 for (i=0;i<p_mc->nI;i++) {
   p_mc->kindI[i]=mc_event->part[i].kind;
   p_mc->momI[i][0]=mc_event->part[i].momentum.t;
   p_mc->momI[i][1]=mc_event->part[i].momentum.space.x;
   p_mc->momI[i][2]=mc_event->part[i].momentum.space.y;
   p_mc->momI[i][3]=mc_event->part[i].momentum.space.z;
 }
 p_mc->nV=mc_event->nfinalVert;
 p_mc->nF=0;
 for (i=0;i<p_mc->nV;i++) {
   mc_vertex=data_getGroup(event,GROUP_MC_EVENT,i+1);
   p_mc->kindV[i]=mc_vertex->kind;
   p_mc->origV[i][0]=mc_vertex->origin.t;
   p_mc->origV[i][1]=mc_vertex->origin.space.x;
   p_mc->origV[i][2]=mc_vertex->origin.space.y;
   p_mc->origV[i][3]=mc_vertex->origin.space.z;	   
   for (j=0;j<mc_vertex->npart;j++) {
     p_mc->kindF[p_mc->nF]=mc_vertex->part[j].kind;
     p_mc->momF[p_mc->nF][0]=mc_vertex->part[j].momentum.t;
     p_mc->momF[p_mc->nF][1]=mc_vertex->part[j].momentum.space.x;
     p_mc->momF[p_mc->nF][2]=mc_vertex->part[j].momentum.space.y;
     p_mc->momF[p_mc->nF][3]=mc_vertex->part[j].momentum.space.z;
     p_mc->vertF[p_mc->nF]=i;
     p_mc->nF++;
   }
 }
 
 return 0;
}

/* * * * * * * * * * * * * * * * * 
* Private method to extract kinematics of given MC photon
* from MC data.
* Returns energy and momentum on success, ow (-1,0) 
*/
static float eventGetPhotonEmc(int *list, float *px, float *py, float *pz) {
 
 int i;
 
 if( (mc_event==NULL) || (p_mc==NULL) ) {
   *px = *py =*pz = 0;	 
   return -1; 
 }
 
 i = *list;
 if ( p_mc->kindF[i] != 1) {
   fprintf(stderr,"ERROR - wrong particle id %d",p_mc->kindF[i]);
   *px = *py =*pz = 0;
   return -1;
 } 
 
 *px = p_mc->momF[i][1];
 *py = p_mc->momF[i][2];
 *pz = p_mc->momF[i][3];
 
 return p_mc->momF[i][0];
 
}

/* * * * * * * * * * * * * * * * * 
* Private method to extract MC photon kinematics
* from the name of the MC itape. Expected format is
*   exxx_txxx_pxxx_.itape
* where xxx should be Energy(e), Theta(t), or Phi(p) value.
*/
static float eventGetPhotonEmc4Name(float *Th, float *Ph) {

 float eng=0;
 float the=0;
 float phi=0;

 char *tok=NULL;
 char *name=NULL;

 if(mcFile==NULL) {
   fprintf(stderr,"ERROR - undefined MC file!\n");
   *Th = *Ph = 0;
   return -1;
 }
 if(name==NULL) {
   name=malloc(1000);
   name[999] = 0;
 }
 strcpy(name,mcFile);

 tok=strtok(name,"/");
 while ( (tok=strtok(NULL,"/")) != NULL ) {
   if (index(tok,'_')) {
     break;
   }
 }

 tok=strtok(tok,"_");
 if (tok[0] == 'e' ) {
   tok++;
   eng = atof(tok);
 }
 else {
   fprintf(stderr,"Error - unknown token %s\n",tok);
   free(name);
   return -1;
 }
 while ( (tok=strtok(NULL,"_")) != NULL ) {
   switch (tok[0]){
     case 't':
       tok++;
       the = atof(tok);
       break;
     case 'p':
       tok++;
       phi = atof(tok);
       break;
     default:
       break;
   }
 }

 *Th = the;
 *Ph = phi;
 free(name);
 name=NULL;

 return eng;
}


/* * * * * * * * * * * * * * * * * 
* Private method to extract angles of isolated showers 
* from input file name.
* The expected form of the name is XXdegYY.itape
* where XX represents theta and YY is slice number).
* The phi=0.2 is default azimuth (11.46 deg).
*/  
static int eventGetPhotonAngle4Name(float *theta, float *phi)  {

 float phi0=THE_PHOTON_PHI;
 
 char *dataFile=NULL;
 char *tok=NULL;
  
 dataFile = malloc(1000);
 if (dataFile==NULL) {
   fprintf(stderr,"ERROR - failed to alocate dataFile\n");
   *theta = *phi = 0;
   return -1;
 }
 dataFile[999]=0;
  
 if (!strcpy(dataFile,inFileName)) {
   fprintf(stderr,"ERROR - copying inFile name %s\n",inFileName);
   *theta = *phi = 0;
   return -1;
 }
 
 tok=strtok(dataFile,"/");
 while ( ((tok=strtok(NULL,"/"))!=NULL) && (strstr(tok,".itape")==NULL) ) { }
 
 tok=strtok(tok,"deg");
 if (tok==NULL) {
   fprintf(stderr,"WARRNING - not an itape with isolated showers!\n");
   *theta = *phi = 0;
   return -1;
 }
 
 *theta = atof(tok)*M_PI/180.;
 *phi = phi0; 
 free(dataFile);
 dataFile=NULL;
 return 0;	
}


/* * * * * * * * * * * * * * * * * 
* Private method to load event characteristic
* (total energy, max energy and channel yields)
* into the lgd_image structire.
* Returns 1 on success, otherwise zero
*/
static int makeLgdImage(void) {
  
 int index=0;
 int channel=0;
 int maxChan=0;
 int row,col;

 float Etot=0;
 float Emax=0;
 float E;
 float x,y;

 memset(lgd_image.yield,0,
        sizeof(float)*LGD_MATRIX_ROWS*LGD_MATRIX_COLS);
 lgd_image.Etot = 0;
 lgd_image.Emax = 0;
 lgd_image.maxChan = 0;

 if((lgdHits==NULL)||(lgdHits->nhits==0)) return 0;

 for (index=0;index<lgdHits->nhits;index++) {
   channel = lgdHits->hit[index].channel;
   E = lgdHits->hit[index].energy;
   Etot += E;
   if (E>Emax) {
     Emax=E;
     maxChan = channel;
   }
   eventGetLgdRC(channel,&row,&col);
   lgd_image.yield[row][col] = E;
 }
 lgd_image.Etot = Etot;
 lgd_image.Emax = Emax;
 lgd_image.maxChan = maxChan;
 eventGetLgdRC(maxChan,&lgd_image.maxRow,&lgd_image.maxCol);

 return 1;  
}

/* * * * * * * * * * * * * * * * *
* Private method to load shower characteristic
* (central channel, total image energy, yields...)
* into the shower_image structure based on hotest LGD block
*/
static int makeShowerImageHot(void) {
   
 int index=0;
 int channel=0;
 int centChan=0;
 int row,col;
 int centRow,centCol;
 int imgRow,imgCol;

 float Etot=0;
 float Emax=0;
 float E;

 int rOffset = SHOWER_IMAGE_OFFSET;
 int cOffset = SHOWER_IMAGE_OFFSET;
 
 memset(lgd_shower_image.yield,0,
        sizeof(float)*SHOWER_IMAGE_DIM*SHOWER_IMAGE_DIM);
 lgd_shower_image.Etot = 0;
 lgd_shower_image.Emax = 0;
 lgd_shower_image.centerChan = 0;
 
 if ((isLgdImg==0) || ((Emax=eventGetLgdEmax(&centChan))==0)) {
   return 0;
 }	 
 
 Emax = eventGetLgdEmax(&centChan);
 eventGetLgdRC(centChan,&centRow,&centCol);

 for(row=(centRow-rOffset);row<=(centRow+rOffset);row++) {
   for(col=(centCol-cOffset);col<=(centCol+cOffset);col++) {
     E=eventGetLgdErc(row,col);
     Etot += E;
     imgRow = row - (centRow-rOffset);
     imgCol = col - (centCol-cOffset);
     lgd_shower_image.yield[imgRow][imgCol] = E;
   }	   
 }
 lgd_shower_image.Etot = Etot;
 lgd_shower_image.Emax = Emax;
 lgd_shower_image.centerChan = centChan;
 eventGetLgdRC(centChan,&lgd_shower_image.centerRow,&lgd_shower_image.centerCol);
 return 1;
} 

/* * * * * * * * * * * * * * * * *
* Private method to load shower characteristic
* (central channel, total image energy, yields...)
* into the shower_image structure based on selected photon 
* from the photon list
*/
static int makeShowerImagePhot(int iphot) {
 
 int index=0;
 int channel=0;
 int centChan=0;
 int row,col;
 int centRow=0,centCol=0;
 int imgRow,imgCol;
 int i;

 float Etot=0;
 float Emax=0;
 float E;
 float x,y;

 double the,phi;
 double xPhot,yPhot;
 double Z0 = LGD_PLANE_DISTANCE; 

 int nbour = SHOWER_CENTER_NEIGHBOUR;
 int rOffset = SHOWER_IMAGE_OFFSET;
 int cOffset = SHOWER_IMAGE_OFFSET;
 
 memset(lgd_shower_image.yield,0,
		        sizeof(float)*SHOWER_IMAGE_DIM*SHOWER_IMAGE_DIM); 
 
 lgd_shower_image.Etot = 0;
 lgd_shower_image.Emax = 0;
 lgd_shower_image.centerChan = 0;

 if((lgdHits==NULL)||(lgdHits->nhits==0)||
    (photons==NULL)||(photons->nPhotons==0)||
    (iPhoton<0)||(isLgdImg==0)) return 0;

/* if there is a photon locate the block closest to it */
 i = iphot;
 phi = atan2(photons->photon[i].momentum.y,
             photons->photon[i].momentum.x);
 the = atan2(sqrt(SQR(photons->photon[i].momentum.x) +
		  SQR(photons->photon[i].momentum.y)),
             photons->photon[i].momentum.z);
 xPhot = Z0 * tan(the) * cos(phi);
 yPhot = Z0 * tan(the) * sin(phi);
 
 
 for (index=0;index<lgdHits->nhits;index++) {
   channel = lgdHits->hit[index].channel;
   E = lgdHits->hit[index].energy;
   eventGetLgdRC(channel,&row,&col);
   eventGetLgdXY(row,col,&x,&y);
   if ( fabs(x-xPhot)<4.0 && fabs(y-yPhot)<4.0 ) {
     centChan = channel;
     Emax = E;
   }
 }
 
/* adjust the central chanel if there is near-by block with higher yield.*/
 eventGetLgdRC(centChan,&centRow,&centCol);
 for (row=(centRow-nbour);row<(centRow+nbour);row++) { 
   for (col=(centCol-nbour);col<(centCol+nbour);col++) {	 
     E = lgd_image.yield[row][col];
     if (E > Emax) {
       Emax = E;
       centChan = col + row*LGD_MATRIX_ROWS;  
     }
   }  
 }
 eventGetLgdRC(centChan,&centRow,&centCol);
 
 for(row=(centRow-rOffset);row<=(centRow+rOffset);row++) {
   for(col=(centCol-cOffset);col<=(centCol+cOffset);col++) {
     E=eventGetLgdErc(row,col);
     Etot += E;
     imgRow = row - (centRow-rOffset);
     imgCol = col - (centCol-cOffset);
     lgd_shower_image.yield[imgRow][imgCol] = E;
   }	   
 }
 lgd_shower_image.centerChan = centChan;
 lgd_shower_image.centerRow = centRow;
 lgd_shower_image.centerCol = centCol; 
 lgd_shower_image.Etot = Etot;
 lgd_shower_image.Emax = Emax;
 
 return 1;
}

/* * * * * * * * * * * * * * * * *
* Private method to return an index of the photon from  
* the Lgd area specified by input angles.
*/
static int eventFindPhoton(float theta, float phi) {

 int i;
 int iphot=0;
 
 float sinx=0;
 float siny=0;
 
 float sinx0 = sin(theta)*cos(phi);
 float siny0 = sin(theta)*sin(phi);
 
 if (photons==NULL) { 
   return -1;
 }		 
  
 for (i=0;i<=photons->nPhotons;i++) {
   sinx = photons->photon[i].momentum.x /
          photons->photon[i].energy;
   siny = photons->photon[i].momentum.y /
	  photons->photon[i].energy;
   if ((fabs(sinx-sinx0)<0.02) && (fabs(siny-siny0)<0.02)) {
     iphot = i;
   }
 } 
 
 return iphot;
}

/* public functions */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* Dummy call to load the module
*/
void eventIO(void) {}		/* c entry point */
void eventio_(void) {}		/* fortran entry */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to obtain and prepare an event 
* so that user can get Lgd data (tot energy, max energy, lgdImage ... )
*
* Note: - keeps open the input file and file with the list of input files 
*         until it hits the end of the last input file. 
*         Call eventclose() if you need to start again!!!
*/
int eventGetNext(void) {

  int ret;
  float theta=0;
  float phi=0;
  
  if (fd == 0) {
    switch (getNextFileIn()) {
      case -1:
        return -1;
      case 0:
        return -1;
    }
    if ( (fd=fileOpen(inFileName,"r")) <= 0) {
      return -1;
    }
  }

  if (event == NULL) {
    fprintf(stderr,"Event buffer size %d\n",BUFSIZE);
    event = malloc(BUFSIZE);
    if (event == NULL) {
      fprintf(stderr,"ERROR - failed to allocating event buffer!\n");
      return -1;
    }
  }

  eventNumber=fileGetEvent(fd,event);  
  switch(eventNumber) {
    case -1:
      return -1;  
    case 0:
      fileClose(fd);
      fd=0;
      return eventGetNext();
  }

  if ( (lgdHits=eventSetup(event)) == NULL ) {
#ifdef VERBOSE
    fprintf(stderr,"Error - failed to setup an event\n");
#endif
    return eventGetNext();
  }
  if ( (photons=eventPhotons(event)) == NULL ) {
#ifdef VERBOSE
    fprintf(stderr,"Error - failed to setup photons\n");
#endif
    return eventGetNext();
  }
  
  if ( (isLgdImg=makeLgdImage()) == 0){
#ifdef VERBOSE
    fprintf(stderr,"Warning - failed to load LgdImage\n");
#endif
    return eventGetNext(); 
  }
  
  if ( (isMC=eventIsMC(event)) == 0 ) {
    if(isFirst) {	  
     fprintf(stderr,"Analysing real data...\n");
     isFirst=0;
    }
#ifdef ISOLATED_SHOWERS    
    if ((eventGetPhotonAngle4Name(&theta,&phi)<0)) {
      fprintf(stderr,"ERROR -  no angle information, exit ... \n");
      return -1;
    }
    if (((iPhoton=eventFindPhoton(theta, phi))<0)){
      fprintf(stderr,"Warning - did not find a photon at (%f,%f)\n",
   		    theta,phi);
      return eventGetNext();
    }
#endif 
  }
  else {  
    if ((isMC==1) && (eventSetP_MC(event) != 0) ) {
      fprintf(stderr,"Error - failed to set p_mc\n");
      return eventGetNext();    
    } 
    if (isFirst) {
      fprintf(stderr,"Analysing MC data...\n");
      isFirst=0;
    }
    iPhoton=0;
  }  
#ifdef ISOLATED_SHOWERS    
    if ((isShowerImg=makeShowerImagePhot(iPhoton))==0) {
#ifdef VERBOSE
      fprintf(stderr,"Warning - failed to load ShowerImagePhot\n");
#endif
      return eventGetNext();
    }
#else
    if ( (isShowerImg=makeShowerImageHot()) == 0 ) {
#ifdef VERBOSE
      fprintf(stderr,"Warning - failed to load ShowerImageHot\n");
#endif
      return eventGetNext();
    }
#endif  
#ifdef VERBOSE
  fprintf(stderr,"Event number: %d\n",eventNumber);
#endif  
  return eventNumber;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* Fortran interfece to get next event
*/
int eventnext_(void) {
 return eventGetNext();
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to close and flush all open files
* 
* Note: - should be called at the and of job!
*/ 
void eventclose_(void) {

  if (inFileName) {
    fprintf(stderr,"Releasing input file name %s ...\n",inFileName);
    free(inFileName);
    inFileName = NULL;
  }	
  if (fList) {
    fprintf(stderr,"Closing list FILE ...\n");
    fclose(fList);
    fList=NULL;
  }
  if (fd) {
    fprintf(stderr,"Closing input file ...\n");	  
    fileClose(fd);
    fd=0;
  }
  return;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to obtain total energy in the LGD shower image
*/
float eventGetEtot() {
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  return lgd_shower_image.Etot;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to obtain total energy in the LGD 
*/
float eventGetLgdEtot() { 
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  return lgd_image.Etot;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to obtain the max energy and corresponding  
* channel of the LGD shower image
*/
float eventGetEmax(int *ch) {
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  *ch = lgd_shower_image.centerChan;
  return lgd_shower_image.Emax;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to obtain the max energy and corresponding  
* channel of the LGD
*/
float eventGetLgdEmax(int *ch) {
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  *ch = lgd_image.maxChan;
  return lgd_image.Emax;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to obtain corresponding row and column 
* from  the given Lgd channel
*/  
void eventGetLgdRC(int ch, int *row, int *col) {
  *row = ch/LGD_MATRIX_ROWS;
  *col = ch - (*row)*LGD_MATRIX_COLS;    
  return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to retrieve XY coordinates, in the Lgd plane,
* of the center of block defined by the given row and column 
*/
void eventGetLgdXY(int row, int col, float *x, float *y) {
   if ( row<0 || row>=LGD_MATRIX_ROWS ||
        col<0 || col>=LGD_MATRIX_COLS ) {
     fprintf(stderr,"WARNING - row/column %d/%d out of LGD,event %d \n",row,col,
		     eventNumber);
     *x = -1000.;
     *y = -1000.;
     return;
   }
   *x = (col - (LGD_MATRIX_COLS/2) + 0.5)*LGD_BLOCK_SIZE;
   *y = (row - (LGD_MATRIX_ROWS/2) + 0.5)*LGD_BLOCK_SIZE;
   return;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to retrieve the energy in a given channel 
* from the lgd_shower_image
*/
float eventGetE(int ch) {
  int row,col;
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  eventGetLgdRC(ch,&row,&col);
  if ( row >= (lgd_shower_image.centerRow-SHOWER_IMAGE_OFFSET) && 
       row <= (lgd_shower_image.centerRow+SHOWER_IMAGE_OFFSET) &&
       col >= (lgd_shower_image.centerCol-SHOWER_IMAGE_OFFSET) && 
       col <= (lgd_shower_image.centerCol+SHOWER_IMAGE_OFFSET) ) {
    row -= lgd_shower_image.centerRow-SHOWER_IMAGE_OFFSET;
    col -= lgd_shower_image.centerCol-SHOWER_IMAGE_OFFSET;
    return lgd_shower_image.yield[row][col];
  }
  return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* User call-able method to retrieve the energy in a given LGD block
*/
float eventGetLgdE(int ch) {
  int row,col;
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  eventGetLgdRC(ch,&row,&col);
  return lgd_image.yield[row][col];
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to retrieve the energy in a given 
* row/col from the lgd_shower_image
*/
float eventGetErc(int row, int col) {
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  if ( row >= (lgd_shower_image.centerRow-SHOWER_IMAGE_OFFSET) && 
       row <= (lgd_shower_image.centerRow+SHOWER_IMAGE_OFFSET) &&
       col >= (lgd_shower_image.centerCol-SHOWER_IMAGE_OFFSET) && 
       col <= (lgd_shower_image.centerCol+SHOWER_IMAGE_OFFSET) ) {
    row -= lgd_shower_image.centerRow-SHOWER_IMAGE_OFFSET;
    col -= lgd_shower_image.centerCol-SHOWER_IMAGE_OFFSET;
    return lgd_shower_image.yield[row][col];
  }
  return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* User call-able method to retrieve energy from the given LGD row and col
*/
float eventGetLgdErc(int row, int col) {
  if (eventNumber <= 0) {
    fprintf(stderr,"ERROR - wrong eventNumber!\n");
    return -1.;
  }
  return lgd_image.yield[row][col];
}

/* * * * * * * * * * * * * * * * * * * * * * * *
* User call-able method to obtain photon four-vector
* found in the slice of around angle defined by THE_PHOTON_PHI
*/
float eventGetPhotonE(int *list,float *Px, float *Py, float *Pz) {

  int i = *list;
  if (i == -1) {
    *Px = *Py = *Pz = 0;
    return -1;
  }  
  *Px = photons->photon[i].momentum.x;
  *Py = photons->photon[i].momentum.y;
  *Pz = photons->photon[i].momentum.z;
  return photons->photon[i].energy;
}



/* * * * * * * * * * * * * * * * * * * * * * * *
* Fortran interface for eventGetPhotonEmc
* At this point, for single MC showers use predefined list-index
* of MC photons (0).
*/
float mcphoton_(float *Px, float *Py, float *Pz) {
	
  int i=0;
  float eng=0;
  float the=0;
  float phi=0;
  
  switch (isMC) {
    case 1:  
      return eventGetPhotonEmc(&i,Px,Py,Pz);
    case 3:
      if ((eng=eventGetPhotonEmc4Name(&the,&phi))<0) {
	*Px = *Py = *Pz = 0;
	return -1;
      }
      the *= (M_PI/180.);
      phi *= (M_PI/180.);
      *Px = eng*sin(the)*cos(phi);
      *Py = eng*sin(the)*sin(phi);
      *Pz = eng*cos(the);
      return eng;
    default:
      *Px = *Py = *Pz = 0;      
      return 0;
  }
}

/* * * * * * * * * * * * * * * * * * * * * * * *
* Fortran interface for eventGetPhotonE
*/
float eventphoton_(float *Px, float *Py, float *Pz) {
  int i = iPhoton;
  return eventGetPhotonE(&i,Px,Py,Pz);
}

/* * * * * * * * * * * * * * * * * 
* Fortran interface for to obtain MC status
*/  
int eventmc_(void) {
 return isMC;
}

/* * * * * * * * * * * * * * * * * * 
* User call-able method to fill a histogram with 
* the image of the LGD hits.
*/
void eventimage_() {
  
  int index;
  int channel;
  int row,col;
  float E;

  hbook2(99,"LGD image ",28,-0.5,27.5,28,-0.5,27.5,0);

  for (row=0;row<LGD_MATRIX_ROWS;row++) {
    for(col=0;col<LGD_MATRIX_COLS;col++) {
      E = eventGetLgdErc(row,col);
      hf2(99,col,row, E); 
    }
  }
  return;
}

/* * * * * * * * * * * * * * * * * * 
* User call-able method to obtain runNumber info
* fortran interface
*/
int  eventrun_() {
  return event->runNo;
}

/* * * * * * * * * * * * * * * * * * 
* User call-able method to call define macros 
* depending on the type of an event. 
* fortran interface
*/
void eventdef_(int *defin) {
 
  define = *defin;	
  return;
}



/* * * * * * * * * * * * * * * * * * 
* User call-able method for testing 
*/
int eventtest_(int *count) {

  int ret=-1;
  int records=0;
  int i=10;
  int index;
  int ch=0;
  int col,row;
  int centCh=0;
  int maxCh=0;
  float Etot=0;
  float Emax=0;
  float EtotLgd=0;
  float EmaxLgd=0; 
  float E,px,py,pz;
  float phi,the;
  float phiRad,theRad;
  int chtest=487;
  float x,y;

  fprintf(stderr," Processing %d events...\n",*count);
  while ( records < *count && (ret=eventGetNext()) > 0 ) {
    records++;
    EtotLgd = eventGetLgdEtot();
    EmaxLgd = eventGetLgdEmax(&ch);
    maxCh = ch;

    Etot = eventGetEtot();
    Emax = eventGetEmax(&ch);
    eventGetLgdRC(ch,&row,&col);
    centCh=ch;

    fprintf(stderr,"LGD: E %f\t Emax %f\t Ch %d\n",EtotLgd,EmaxLgd,maxCh);
    fprintf(stderr,"IMG: E %f\t Emax %f\t Ch %d\n",Etot,Emax,centCh);
    hf1(20,EtotLgd,1.);
    hf2(21,col,row,EmaxLgd);
    hf1(22,Etot,1.);
    hf2(23,col,row,Emax);
    for (index=0;index<ndim;index++) {
      eventGetLgdRC(index,&row,&col);
      hf2(200,col,row,eventGetErc(row,col));
    }

#ifdef ISOLATED_SHOWERS
    if ((ret=eventGetPhotonAngle4Name(&theRad,&phiRad)) < 0) {
      fprintf(stderr,"Warning - can't get angle from input file name %sn",
	      inFileName);
    };

    if ((iPhoton=eventFindPhoton(theRad,phiRad)) <0) {
       fprintf(stderr,"Error - negative photon index %d\n",iPhoton);
    }
#endif
    
    E = eventGetPhotonE(&iPhoton,&px,&py,&pz);
    phi = atan2(py,px);
    the = atan2(sqrt(SQR(px)+SQR(py)),pz);
    fprintf(stderr,"Photon %d\t E %f\t theta %f\t phi %f\n",
                 iPhoton,E,the,phi);

    if(isMC) {
      E = mcphoton_(&px,&py,&pz);
      phi = atan2(py,px);
      the = atan2(sqrt(SQR(px)+SQR(py)),pz);
      fprintf(stderr,"MC photon: E %f\t theta %f\t phi %f\n",
                                     E,the,phi);
    }

    chtest=ch+1;
    eventGetLgdRC(chtest,&row,&col);
    eventGetLgdXY(row,col,&x,&y);
    fprintf(stderr,"Ch %d, Row %d, Col %d, X %f, Y %f, Yield %f \n",
                  chtest,row,col,x,y,eventGetErc(row,col));

    if (records/i*i == records) { 
      fprintf(stderr,"events %d, event %d\n",records,ret);
    }
  }
  return 0;
}

#include "eventIOfriends.c"

