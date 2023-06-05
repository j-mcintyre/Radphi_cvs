/*
 * March 10, 2004
 *   Modifications to data selection cuts:
 *    - eliminated MAX_CPVIN_TIME_HITS cut as incompatible with tagging
 *    - introduced MAKE_ALL to make all groups
 *    - support for new GROUP_RECOIL
 *    - support for new GROUP_TAGGER_CLUSTERS
 *    - support for new GROUP_TAGGER_PHOTONS
 *    - new cut on neutral BGV energy
 */
/* 
 * February 25, 2003
 *   Introduced data selection cuts: 
 *    - number of LGD photons
 *    - minimum cluster separation
 *    - total energy in the LGD
 *    - number of BSD pixels
 *    - CPV hits within the time window from given offset 
 *    - mesons selection
 *    - minimum cluster energy
 *   For each selection criterion choose corresponding macro and value.
 *   The total LGD and cluster energy belong to continuum, thus the 
 *   default cut values can be superseded from the command line. The 
 *   discrete cuts are currently hard-wired. Nevertheless, in order for 
 *   selection cut to have an effect corresponding macro should be enabled.
 * - MK
*/    
/*
 * February 14, 2003
 *   Added ability to save all events that do not hit the target spots
 *   to another output file called 99deg.itape, selected by the macro
 *   SAVEALL.
 * - Richard Jones
 */

//#define CLUSTER_CLEANUP		      

/* 
 *  Select at least one method below, in order to produce output. 
 */
//#define SELECT_ISOLATED_SHOWERS true
//#define SAVEALL true
#define USE_EVENT_CUTS true

/* 
 * Select at least one event selection criteria below,
 * if you have enabled USE_EVENT_CUTS above.
 */
#define NUM_LGD_PHOTONS 7
#define MIN_TOTAL_LGD_ENERGY 3.0	// GeV
#define MIN_CLUSTER_ENERGY 0.05		// GeV
//#define MIN_CLUSTER_SEPARATION 15	// cm
//#define NUM_BSD_PIXELS 1
//#define USE_MESON 66 			// meson type 
//#define NUM_MESONS 1                  
#define NUM_RECOILS 1
//#define NUM_NEUTRAL_BGV 0
//#define NUM_BGV_HITS 1
//#define NUM_CPV_COINS 0
/* In order to use FULLY_CONTAINED_FORWARD one has to enable
 * at least neutral bgv and recoils out of the previous 3 filters */
//#define FULLY_CONTAINED_FORWARD 1

#define MAXIMUM_RECOIL_TIME 10.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <math.h>

#include <dataIO.h>
#include <param.h>
#include <triggerType.h>
#include <makeHits.h>
#include <lgdGeom.h>
#include <libBGV.h>
#include <umap.h>

#define CPV_RECOIL_WINDOW_LOW -3.0	// ns
#define CPV_RECOIL_WINDOW_HIGH +3.0	// ns
#define LGD_PLANE_DISTANCE 120

#define SQR(X) ((X)*(X))

FILE *fpout[]={0,0,0,0,0,0,0};

double bytes=0;
int records=0;
int nData=0;
int nTrigData=0;
int nTrigMC=0;
int nTrigLGDMon=0;
int nTrigBGVMon=0;
int nTrigBT=0;
int nTrigPed=0;
int nTrigUnknown=0;
int nSync=0;
int nPreStart=0;
int nGo=0;
int nPause=0;
int nEnd=0;
int nUnknown=0;

float EtotThreshold=0.;
float EclustThreshold=0.;

void Usage(void);
static float getMinClustSep(photons_t *photons);
static float getMinClustEnergy(photons_t *photons);
static float getEtotLGD(lgd_hits_t *lgd_hits);
static int getMesonId(mesons_t *mesons);
static int getNbgvNeutral(bgv_hits_t *bgv_hits,bsd_pixels_t *bsd_pixs);
static int getCPvetos(cpv_hits_t *cpv_hits, recoils_t *recoils,
                      time_list_t *time_list);
static int getNrecoils(recoils_t *recoils);

int psfReduce (char *inspec)
{
  int *raw;
  int ret=0;
  int srun[]={0,0};
  static itape_header_t *event=NULL;
  int runNo=999999;
  char infile[1000];
  char *fin;
  int fdin;
  int i;
  int nphotons_in_lgd;

  mesons_t *mesons=NULL;
  photons_t *photons=NULL;
  lgd_hits_t *lgd_hits=NULL;
  bsd_pixels_t *bsd_pixs=NULL;
  time_list_t *time_list=NULL;
  bgv_hits_t *bgv_hits=NULL;
  cpv_hits_t *cpv_hits=NULL;
  recoils_t *recoils=NULL;
   
  
/* Malloc the event buffer */
  if (event == NULL) {
    event = malloc(BUFSIZE);
    if (event == NULL) {
      fprintf(stderr,"ERROR - failed to allocating event buffer!\n");
      return -1;
    }
  }

/* open the event file */
  strncpy(infile,inspec,1000);
  if ((fin=(char *)index(infile,' ')) != NULL) {
    fin[0]=0;
  }
  if ((fdin=evt_open(infile,"r")) == 0) {
    fprintf(stderr,"ERROR - failed to open input event file \"%s\"\n",infile);
    return -1;
  }

#if defined SELECT_ISOLATED_SHOWERS
  for (i=0; i<5; i++) {
    if (fpout[i] == 0) {
      char filename[99];
      sprintf(filename,"%ddeg.itape",5*(i+1));
      if ((fpout[i] = fopen(filename,"w")) == 0) {
        fprintf(stderr,"ERROR - failed to open output file %s!\n",filename);
        return -1;
      }
      else {
        tape_setStream(fpout[i]);
      }
    }
    data_setSieve(fpout[i],fdin,infile);
  }
#endif 
#if defined SAVEALL
  if (fpout[5] == 0) {
    char filename[99];
    sprintf(filename,"99deg.itape");
    if ((fpout[5] = fopen(filename,"w")) == 0) {
      fprintf(stderr,"ERROR - failed to open output file!\n");
      return -1;
    }
    else {
      tape_setStream(fpout[5]);
    }
  } 
  data_setSieve(fpout[5],fdin,infile);
#endif
#if defined USE_EVENT_CUTS
  if (fpout[6] == 0) {
    char filename[99];   
    sprintf(filename,"sample.itape");
    if ((fpout[6] = fopen(filename,"w")) == 0) {
      fprintf(stderr,"ERROR - failed to open output file!\n");
      return -1; 
    }
    else {
      tape_setStream(fpout[6]);
    } 
  }	 
  data_setSieve(fpout[6],fdin,infile); 
#endif  
  
/* loop over input events */
  while  ((ret=evt_data_read(fdin,event,BUFSIZE)) == DATAIO_OK) {
    records++;
    bytes += event->length + 4;
    if (event->eventType == EV_DATA) {
      if (event->runNo != srun[0]) {
        if (setup_makeHits(event)) {
           fprintf(stderr,"ERROR - setup_makeHits failed!\n");
           return -1;
        }
        srun[0] = event->runNo;
        srun[1] = 0;
      }
      /* setup stuff */
      if (event->trigger == TRIG_DATA) {
        unpackEvent(event,BUFSIZE);
        makeHits(event);
      }
      data_clean(event);
      /* event selection */
#if defined SELECT_ISOLATED_SHOWERS	
      if (
#if defined CLUSTER_CLEANUP		      
	  (lgd_cluster_cleanup(event) == 0) &&
#endif		      
	 (photons=data_getGroup(event,GROUP_PHOTONS,0)) != NULL) {
        int i;
        double phi0=0.2;
        double sinx0 = cos(phi0)*sin(5 * M_PI/180);
        double siny0 = tan(phi0)*sinx0;
        double sinx1 = cos(phi0)*sin(10 * M_PI/180);
        double siny1 = tan(phi0)*sinx1;
        double sinx2 = cos(phi0)*sin(15 * M_PI/180);
        double siny2 = tan(phi0)*sinx2;
        double sinx3 = cos(phi0)*sin(20 * M_PI/180);
        double siny3 = tan(phi0)*sinx3;
        double sinx4 = cos(phi0)*sin(25 * M_PI/180);
        double siny4 = tan(phi0)*sinx4;
        for (i=0; i<photons->nPhotons; i++) {
          double sinx = photons->photon[i].momentum.x /
                        photons->photon[i].energy;
          double siny = photons->photon[i].momentum.y /
                        photons->photon[i].energy;
          double sinxx, sinyy;
          int ii;
          for (ii=0; ii<photons->nPhotons; ii++) {
            if (ii == i) { continue; }
            sinxx = photons->photon[ii].momentum.x /
                    photons->photon[ii].energy;
            sinyy = photons->photon[ii].momentum.y /
                    photons->photon[ii].energy;
            if ((fabs(sinx-sinxx)<0.20) && (fabs(siny-sinyy)<0.20)) {
              break;
            }
          }
          if (ii < photons->nPhotons) {
            continue;
          }
          if ((fabs(sinx-sinx0)<0.01) && (fabs(siny-siny0)<0.01)) {
            data_write(fileno(fpout[0]),event);
          }
          else if ((fabs(sinx-sinx1)<0.01) && (fabs(siny-siny1)<0.01)) {
            data_write(fileno(fpout[1]),event);
          }
          else if ((fabs(sinx-sinx2)<0.01) && (fabs(siny-siny2)<0.01)) {
            data_write(fileno(fpout[2]),event);
          }
          else if ((fabs(sinx-sinx3)<0.01) && (fabs(siny-siny3)<0.01)) {
            data_write(fileno(fpout[3]),event);
          }
          else if ((fabs(sinx-sinx4)<0.01) && (fabs(siny-siny4)<0.01)) {
            data_write(fileno(fpout[4]),event);
          }
        }
      }
#endif
#if defined SAVEALL
      data_write(fileno(fpout[5]),event);
#endif

      nphotons_in_lgd=0;
      if ((photons=data_getGroup(event,GROUP_PHOTONS,0)) != NULL) {
         for (i=0; i < photons->nPhotons; i++) {
            if (photons->photon[i].cluster >= 0) {
               ++nphotons_in_lgd;
            }
         }
      }

#if defined USE_EVENT_CUTS
      if (
#  if defined CLUSTER_CLEANUP		      
		   (lgd_cluster_cleanup(event) == 0) &&
#  endif		      
          ((lgd_hits=data_getGroup(event,GROUP_LGD_HITS,0)) != NULL) &&
	  ((photons=data_getGroup(event,GROUP_PHOTONS,0)) != NULL) && 
	  getEtotLGD(lgd_hits) >= EtotThreshold && 
	  getMinClustEnergy(photons) >= EclustThreshold )
# if defined NUM_LGD_PHOTONS
      if (((photons=data_getGroup(event,GROUP_PHOTONS,0)) != NULL) && 
          nphotons_in_lgd == NUM_LGD_PHOTONS )
# endif
# if defined MIN_CLUSTER_SEPARATION
      if (((photons=data_getGroup(event,GROUP_PHOTONS,0)) != NULL) &&  
	  getMinClustSep(photons) >= MIN_CLUSTER_SEPARATION)
# endif
# if defined NUM_BSD_PIXELS 
      if (((bsd_pixs=data_getGroup(event,GROUP_BSD_PIXELS,0)) != NULL) &&
           bsd_pixs->npixels == NUM_BSD_PIXELS)
# endif
# if defined NUM_CPV_COINS
      if (((cpv_hits=data_getGroup(event,GROUP_CPV_HITS,0)) != NULL) &&
          ((time_list=data_getGroup(event,GROUP_TIME_LIST,0)) != NULL) &&
          ((recoils=data_getGroup(event,GROUP_RECOIL,0)) != NULL) &&
	   (recoils->nrecoils > 0) &&
           getCPvetos(cpv_hits,recoils,time_list) == NUM_CPV_COINS)
# endif
# if defined USE_MESON
      if (((mesons=data_getGroup(event,GROUP_MESONS,0)) != NULL) &&
           getMesonId(mesons) == USE_MESON)
# endif
# if defined FULLY_CONTAINED_FORWARD && NUM_NEUTRAL_BGV && NUM_RECOILS
      if (((bgv_hits=data_getGroup(event,GROUP_BGV_HITS,0)) != NULL) 
	  && ((bsd_pixs=data_getGroup(event,GROUP_BSD_PIXELS,0)) != NULL) 
	  && ((recoils=data_getGroup(event,GROUP_RECOIL,0)) != NULL) 
	  && (getNbgvNeutral(bgv_hits,bsd_pixs) == NUM_NEUTRAL_BGV)
	  && (getNrecoils(recoils) == NUM_RECOILS)
#   if defined NUM_BGV_HITS
	  &&  (bgv_hits->nhits <= NUM_BGV_HITS)
#   endif	   
       	               ) 
# else
#   if defined NUM_BGV_HITS
      if (((bgv_hits=data_getGroup(event,GROUP_BGV_HITS,0)) != NULL) &&
	   bgv_hits->nhits <= NUM_BGV_HITS)	      
#   endif
#   if defined NUM_NEUTRAL_BGV
      if (((bgv_hits=data_getGroup(event,GROUP_BGV_HITS,0)) != NULL) &&
	  ((bsd_pixs=data_getGroup(event,GROUP_BSD_PIXELS,0)) != NULL) &&
	    getNbgvNeutral(bgv_hits,bsd_pixs) == NUM_NEUTRAL_BGV)	      
#   endif
#   if defined NUM_RECOILS
      if (((recoils=data_getGroup(event,GROUP_RECOIL,0)) != NULL) &&
	  getNrecoils(recoils) == NUM_RECOILS)
#   endif	      
# endif	      
        data_write(fileno(fpout[6]),event);
#endif
      ++srun[1];
    }
    switch(event->eventType) {
    case EV_DATA:
      switch(event->trigger) {
      case TRIG_DATA:
        nTrigData++;
	break;
      case TRIG_MC:
        nTrigMC++;
	break;
      case TRIG_LGDMON:
	nTrigLGDMon++;
	break;
      case TRIG_PED:
        nTrigPed++;
        break;
      case TRIG_BT:
        nTrigBT++;
        break;
      case TRIG_BGVMON:
        nTrigBGVMon++;
        break;
      default:
        nTrigUnknown++;
        break;
      }
      nData++;
      break;
    case EV_SYNC:
      nSync++;
      break;
    case EV_PRESTART:
      nPreStart++;
      break;
    case EV_GO:
      nGo++;
      break;
    case EV_PAUSE:
      nPause++;
      break;
    case EV_END:
      nEnd++;
      break;
    case EV_UNKNOWN:
      nUnknown++;
      break;
    default:
      fprintf(stderr,"EVent type : %d\n",event->eventType);
    }	  
    {
      int i;
      for (i=10000;i<records;i*=10) {}
      i/=10;
      if (records/i*i == records) {
        printf("%d\n",records);
      }
    }
  }
  if (ret == DATAIO_BADCRC || ret == DATAIO_ERROR) {
     fprintf(stderr,"error occurred in event %d after %f bytes read.\n",
             records,bytes);
  }
  {
     int i;
     for (i=0;i<7;i++) {
        if (fpout[i]) {
           data_flush(fileno(fpout[i]));
        }
     }
  }
  evt_close(fdin);
  return (ret < 0)? ret:0;
}

/*
 * keep this around -- it limits the size of core dumps!
    struct rlimit rlim;
    getrlimit(RLIMIT_CORE,&rlim);
    rlim.rlim_cur = 1000;
    setrlimit(RLIMIT_CORE,&rlim);
 */

int main(int argc, char *argv[])
{
  int arg;
  int ret=0;
#if defined MIN_TOTAL_LGD_ENERGY
  EtotThreshold = MIN_TOTAL_LGD_ENERGY;
#endif
#if defined MIN_CLUSTER_ENERGY
  EclustThreshold = MIN_CLUSTER_ENERGY;
#endif
  for (arg=1; arg < argc; arg++) {
    if (*argv[arg] == '-') {
      char *argp = argv[arg]+1;
      switch (*argp) {
        case 'E':
	  EtotThreshold = atof(argv[++arg]);
	  break;
	case 'e':
	  EclustThreshold = atof(argv[++arg]);
	  break;
	case 'h':
	  Usage();
          exit(1);
	default:  
	  fprintf(stderr,"Unknown argument: %s\n",argv[arg]);
	  Usage();
	  exit(1);
      }
    }
    else {
      ret = psfReduce(argv[arg]);
      if (ret != 0) {
        perror("quitting");
        break;
      }
    }	    
  }

  fprintf(stderr,"Processed %d records:\n",records);
  fprintf(stderr,"\t%d data",nData);
  fprintf(stderr,"\t%d physics",nTrigData);
  fprintf(stderr,"\t%d simulation",nTrigMC);
  fprintf(stderr,"\t%d LGD monitor",nTrigLGDMon);
  fprintf(stderr,"\t%d BGV monitor\n",nTrigBGVMon);
  fprintf(stderr,"\t\t\t%d pedestal",nTrigPed);
  fprintf(stderr,"\t%d basetest",nTrigBT);
  fprintf(stderr,"\t%d unknown\n",nTrigUnknown);
  fprintf(stderr,"\t%d sync",nSync);
  fprintf(stderr,"\t%d prestart",nPreStart);
  fprintf(stderr,"\t%d go",nGo);
  fprintf(stderr,"\t%d pause",nPause);
  fprintf(stderr,"\t%d end",nEnd);
  fprintf(stderr,"\t%d unknown\n",nUnknown);
  fprintf(stderr,"%f total bytes read\n",bytes);
  exit(ret);
}


/*
 * Private method to obtain minimum distance between any 2 clusters.
 */ 
static float getMinClustSep(photons_t *photons) {
 
 int i;
 int ii;
 float mdist=1000;
 float dist;
 float Z0 = LGD_PLANE_DISTANCE;
 double xi, yi;
 double xii, yii;
 
	 
 if ((photons == NULL) || (photons->nPhotons < 2)) return -1;	

 for (i=0; i<photons->nPhotons; i++) {
   xi = photons->photon[i].momentum.x/
	photons->photon[i].momentum.z;   
   yi = photons->photon[i].momentum.y/
	photons->photon[i].momentum.z;   
   for (ii=i+1; ii<photons->nPhotons; ii++) {
     xii = photons->photon[ii].momentum.x/
	   photons->photon[ii].momentum.z;   
     yii = photons->photon[ii].momentum.y/
	   photons->photon[ii].momentum.z;   
     dist = Z0*sqrt(SQR(xi-xii) + SQR(yi-yii));
     if (dist < mdist) { mdist=dist; }
   }
   
 }
 return mdist;
}

/*
 * Private method to obtain minimum cluster energy.
 */ 
static float getMinClustEnergy(photons_t *photons) {
 int i;
 float Emin=100.;

 if (photons == NULL) return -1;

 for (i=0; i<photons->nPhotons; i++) {
   float E = photons->photon[i].energy;	 
   if (E < Emin) {Emin=E ;}
 }
 return Emin;
}

/*
 * Private Method to obtain total energy of hits in the LGD. 
 */
static float getEtotLGD(lgd_hits_t *lgd_hits){
	
 int i;
 float E=0;
 if (lgd_hits == NULL) return 0;  
  
 for (i=0;i<lgd_hits->nhits;i++) {
   E += lgd_hits->hit[i].energy;
 }	 
 return E;
}


static int getMesonId(mesons_t *mesons) {

 int i;
 int id;
 int nid=0;


 if (mesons == NULL) return -1;	
#if defined USE_MESON 
 for (i=0; i<mesons->nMesons; i++) {
   if ((id=mesons->meson[i].type) == USE_MESON) {
     nid++;
   }
 }
# if defined  NUM_MESONS 
 if (nid == NUM_MESONS) {return USE_MESON;}
 return -1;
# else
 if (nid) return USE_MESON; 
# endif
#endif 
 return -1; 
}

static int getNbgvNeutral(bgv_hits_t *bgv_hits,bsd_pixels_t *bsd_pixs) {
	
   int neutral_hits=0;
   int i,j,s1, s2;
   for (i=0; i<bgv_hits->nhits; i++) {
     ++neutral_hits;
     for (j=0; j<bsd_pixs->npixels; j++) {
         s1=(bsd_pixs->pixel[j].straight-22)%24;
         s2=(bsd_pixs->pixel[j].straight-23)%24;
	if ( bgv_hits->hit[i].channel == s1 ||
             bgv_hits->hit[i].channel == s2 ) {
	   --neutral_hits; 
	   break;
	}
     }  
   }
   return neutral_hits;
}

static int getCPvetos(cpv_hits_t *cpv_hits, 
                      recoils_t *recoils,
                      time_list_t *time_list) {

   int coins=0;
   float t0;
   int ic;
   if (recoils->nrecoils == 0) {
      return 0;
   }
   else {
      t0 = recoils->best_t0;
   }
   for (ic=0; ic < cpv_hits->nhits; ic++) {
      int it;
      for (it=0; it < cpv_hits->hit[ic].times; it++) {
         int i=it+cpv_hits->hit[ic].tindex;
         float dt=time_list->element[i].le-t0;
         if ((dt > CPV_RECOIL_WINDOW_LOW) && (dt < CPV_RECOIL_WINDOW_HIGH)) {
            coins++;
         }
      }
   }
   return coins;
}

static int getNrecoils(recoils_t *recoils)
{
   int recount=0;
   int irec;
   float t0;
   if (recoils->nrecoils == 0) {
      return 0;
   }
   else {
      t0 = recoils->best_t0;
   }
   for (irec=0; irec < recoils->nrecoils; ++irec) {
      if (recoils->recoil[irec].time-t0 < MAXIMUM_RECOIL_TIME) {
         ++recount;
      }
   }
   return recount;
}

void Usage(void) {
 
 fprintf(stderr,"Usage:\n");
 fprintf(stderr,"  psfReduce <options> file1 file2 ...\n");
 fprintf(stderr,"Options are:\n");
 fprintf(stderr,"\t-E Etot\t\tTotal LGD energy threshold \n");
 fprintf(stderr,"\t-e Eclust\t\tCluster energy threshold \n");
 fprintf(stderr,"\t-h\t\tPrint this message\n");
}
