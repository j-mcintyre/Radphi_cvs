#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <dataIO.h>
#include <param.h>
#include <unpackData.h>
#include <triggerType.h>
#include <makeHits.h>
#include <lgdGeom.h>
#include <libBGV.h>
#include <umap.h>
#include <cernlib.h>

#include "ntuple.h"
#include "hituple.h"

//#define HISTOGRAM_ONLY 1
//#define CLUSTER_CLEANUP 1
//#define FULLY_CONTAINED_FORWARD 1

#define DEFAULT_INIT_PATTERN "arbxvcutlpmoe"

#define sqr(x) ((x)*(x))
#define OUTPUT_FILE_FULL 95
#define OUTPUT_FILE_OVERRUN 96

void merge_lowmass_cluster_pairs(photons_t *photons);

extern struct {
  int iq[100];
} quest_;

int fd=0;

int ntuple_ (char *inspec, int *ntid, int *count, int *start, int len)
{
  int ret=0;
  int srun[]={0,0};
  static itape_header_t *event=NULL;
  int runNo=999999;
  char infile[1000];
  char *fin;
  int i,j;
  int cleanlgd;
  
  static nt_head_t head;
  static nt_head_t *p_head=&head;
  static nt_raw_event_t *p_raw;
  static nt_bsd_hits_t *p_bsd=NULL;
  static nt_bsd_pixs_t *p_pix=NULL;
  static nt_bgv_hits_t *p_bgv=NULL;
  static nt_cpv_hits_t *p_cpv=NULL;
  static nt_upv_hits_t *p_upv=NULL;
  static nt_tag_hits_t *p_tag=NULL;
  static nt_lgd_hits_t *p_lgd=NULL;
  static nt_time_list_t *p_tlist=NULL;
  static nt_photons_t *p_clust=NULL;
  static nt_bclusters_t *p_bclust=NULL;
  static nt_mesons_t *p_mes=NULL;
  static nt_recoils_t *p_rec=NULL;
  static nt_tag_coin_t *p_coin=NULL;
  static nt_tag_align_t *p_talign=NULL;

  int *raw;

  bsd_hits_t *bsd_hits=NULL;
  bsd_pixels_t *bsd_pixs=NULL;
  bgv_hits_t *bgv_hits=NULL;
  cpv_hits_t *cpv_hits=NULL;
  upv_hits_t *upv_hits=NULL;
  lgd_hits_t *lgd_hits=NULL;
  tagger_hits_t *tagger_hits=NULL;
  time_list_t *time_list=NULL;
  lgd_clusters_t *lgd_clusters=NULL;
  bgv_clusters_t *bgv_clusters=NULL;
  photons_t *photons=NULL;
  mesons_t *mesons=NULL;
  recoils_t *recoils=NULL;
  tagger_photons_t *tagging=NULL;

  reg_values_t *mam_regs;
  adc_values_t *bsd_adcs;
  tdc_values_t *bsd_tdcs;
  adc_values_t *bgv_adcs;
  tdc_values_t *bgv_tdcs;
  adc_values_t *cpv_adcs;
  tdc_values_t *cpv_tdcs;
  adc_values_t *upv_adcs;
  tdc_values_t *upv_tdcs;
  tdc_values_t *tagger_tdcs;
  adc_values_t *lgd_adcs;

  int bytes=0;
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

/* Malloc the event buffer */
  if (event == NULL) {
    event = malloc(BUFSIZE);
    if (event == NULL) {
      fprintf(stderr,"Error - failed to allocating event buffer!\n");
      return -1;
    }
  }

  if (strncmp(inspec,"init",4) == 0) {
    char *tok = malloc(100);
    struct rlimit rlim;
    getrlimit(RLIMIT_CORE,&rlim);
    rlim.rlim_cur = 1000;
    setrlimit(RLIMIT_CORE,&rlim);
    if (p_raw)   free(p_raw), p_raw = NULL;
    if (p_bsd)   free(p_bsd), p_bsd = NULL;
    if (p_pix)   free(p_pix), p_pix = NULL;
    if (p_bgv)   free(p_bgv), p_bgv = NULL;
    if (p_cpv)   free(p_cpv), p_cpv = NULL;
    if (p_upv)   free(p_upv), p_upv = NULL;
    if (p_tag)   free(p_tag), p_tag = NULL;
    if (p_lgd)   free(p_lgd), p_lgd = NULL;
    if (p_tlist) free(p_tlist), p_tlist = NULL;
    if (p_clust) free(p_clust), p_clust = NULL;
    if (p_bclust) free(p_bclust), p_bclust = NULL;
    if (p_mes)   free(p_mes), p_mes = NULL;
    if (p_rec)   free(p_rec), p_rec = NULL;
    if (p_coin)  free(p_coin), p_coin = NULL;
    if (p_talign) free(p_talign), p_talign = NULL;

/* decide which of the ntuple blocks should be included by
 * looking at the string of letters following the init key.
 */
    strncpy(tok,inspec,100);
    tok[99] = 0;
    tok = strtok(tok," ");
    if ((tok = strtok(NULL," ")) == NULL) {
      tok = malloc(80);
      sprintf(tok,DEFAULT_INIT_PATTERN);
    }
    if (index(tok,'r') != NULL) {
      p_raw = malloc(sizeof(nt_raw_event_t));
    }
    if (index(tok,'b') != NULL) {
      p_bsd = malloc(sizeof(nt_bsd_hits_t));
    }
    if (index(tok,'x') != NULL) {
      p_pix = malloc(sizeof(nt_bsd_pixs_t));
    }
    if (index(tok,'v') != NULL) {
      p_bgv = malloc(sizeof(nt_bgv_hits_t));
    }
    if (index(tok,'c') != NULL) {
      p_cpv = malloc(sizeof(nt_cpv_hits_t));
    }
    if (index(tok,'u') != NULL) {
      p_upv = malloc(sizeof(nt_upv_hits_t));
    }
    if (index(tok,'t') != NULL) {
      p_tag = malloc(sizeof(nt_tag_hits_t));
    }
    if (index(tok,'l') != NULL) {
      p_lgd = malloc(sizeof(nt_lgd_hits_t));
    }
    if (p_bsd || p_bgv || p_cpv || p_upv || p_tag) {
      p_tlist = malloc(sizeof(nt_time_list_t));
    }
    if (index(tok,'a') != NULL) {
      p_bclust = malloc(sizeof(nt_bclusters_t));
    }
    if (index(tok,'p') != NULL) {
      p_clust = malloc(sizeof(nt_photons_t));
    }
    if (index(tok,'m') != NULL) {
      p_mes = malloc(sizeof(nt_mesons_t));
    }
    if (index(tok,'e') != NULL) {
      p_rec = malloc(sizeof(nt_recoils_t));
    }
    if (index(tok,'o') != NULL) {
      p_coin = malloc(sizeof(nt_tag_coin_t));
    }
    if (index(tok,'i') != NULL) {
      p_talign = malloc(sizeof(nt_tag_align_t));
    }

    return 0;
  }

/* open the event file */
  strncpy(infile,inspec,1000);
  if ((fin=(char *)index(infile,' ')) != NULL)
    fin[0] = 0;
  if ((fd == 0) && ((fd=evt_open(infile,"r")) == 0)) {
    fprintf(stderr,"Error - failed to open input event file \"%s\"\n",
            infile);
    return -1;
  }

  if (*start > 0) {
    printf("Skipping %d events...\n",*start);
    while (*start) {
      ret = evt_data_read(fd,event,BUFSIZE);
      switch (ret) {
      case DATAIO_OK:
        bytes += event->length + 4;
        records++;
        (*start)--;
        break;
      case DATAIO_BADCRC:
      case DATAIO_ERROR:
        fprintf(stderr,"Error occurred in event %d after %d bytes read.\n",
                records,bytes);
      default:
        evt_close(fd);
        fd = 0;
        return ret;
      }
    }
  }

  while (*count) {
    ret = evt_data_read(fd,event,BUFSIZE);
    if (ret != DATAIO_OK)
      break;
    bytes += event->length + 4;
    records++;
    (*count)--;

    switch (event->eventType) {
    case EV_DATA:
      if (event->runNo != srun[0]) {
        if (setup_makeHits(event)) {
          fprintf(stderr,
                  "Error in ntuple: setup_makeHits failed, quitting!\n");
          return -1;
        }
        else {
          char ccFile[] = "initial_lgd_cc_values.dat";
          FILE *fpin = fopen(ccFile,"r");
          if (fpin) {
            int col,row,chan;
            float gain,gc;
            int index;
            for (index=0; index < lgd_nChannels; index++) {
              fscanf(fpin,"%d %d %d %f %f",
                     &chan,&col,&row,&gain,&gc);
              if (chan != index){
                fprintf(stderr,"cc input error: ");
                fprintf(stderr,"channel %d out of order",chan);
                fprintf(stderr," in input file %s\n",ccFile);
                exit(1);
              }
              else {
                lgd_cc[index]*=gain;
              }
            }
            fclose(fpin);
          }
        }
        srun[0] = event->runNo;
        srun[1] = 0;
      }

      unpackEvent(event,BUFSIZE);

      p_head->runNo = event->runNo;
      p_head->eventNo = event->eventNo;
      p_head->isMC = (event->trigger == TRIG_MC);

      if (p_raw) {
        p_raw->evType = event->trigger;
  
        if ((mam_regs=data_getGroup(event,GROUP_MAM_REGS,0)) != NULL) {
           p_raw->meMAM = mam_regs->reg[0].contents;
           p_raw->cMAM = mam_regs->reg[4].contents;
           p_raw->eMAM = mam_regs->reg[7].contents;
           p_raw->m2MAM = mam_regs->reg[8].contents;
        }
  
        p_raw->nBSD = 0;
        if ((bsd_tdcs=data_getGroup(event,GROUP_BSD_TDCS_LONG,0)) != NULL) {
          for (i=0;i<bsd_tdcs->ntdc;i++) {
            if ((bsd_tdcs->tdc[i].le >= TMIN_BSD) &&
                (bsd_tdcs->tdc[i].le <= TMAX_BSD)) {
              if (p_raw->nBSD > NT_MAX_BSD_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in bsd, maximum allowed was %d\n",
                          p_raw->nBSD, NT_MAX_BSD_HIT);
                break;
              }
              else {
                p_raw->iBSD[p_raw->nBSD] = bsd_tdcs->tdc[i].channel;
                p_raw->tBSD[p_raw->nBSD] = bsd_tdcs->tdc[i].le;
                p_raw->aBSD[p_raw->nBSD] = 0;
                p_raw->nBSD++;
              }
            }
          }
        }
  /* For counters like the BSD that are read out with multihit TDCs
   * I have assigned a special meaning to the tBSD() array.  Each
   * TDC hit creates its own entry in the xBSD() arrays, with the ADC
   * (if any) duplicated for the mulitple TDC values.  This is quite
   * convenient for studying time spectra, but messes up the statistics
   * of ADC spectra.  For this reason, I add an extra entry at the end
   * of the xBSD() arrays for each unique aBSD, for which the tBSD
   * contains the count of TDC values for that channel instead of a
   * TDC value.  For studies of timing spectra, the low integers in
   * tBSD should be avoided.                                RTJ  July 1999
   */
        if ((bsd_adcs=data_getGroup(event,GROUP_BSD_ADCS,0)) != NULL) {
          for (i=0;i<bsd_adcs->nadc;i++) {
            if (bsd_adcs->adc[i].value >= AMIN_BSD) {
              int ntdcs=0;
              for (j=0;j<p_raw->nBSD;j++) {
                if (p_raw->iBSD[j] == bsd_adcs->adc[i].channel) {
                  p_raw->aBSD[j] = bsd_adcs->adc[i].value;
                  ntdcs++;
                }
              }
              if (p_raw->nBSD > NT_MAX_BSD_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in bsd, maximum allowed was %d\n",
                          p_raw->nBSD, NT_MAX_BSD_HIT);
                break;
              }
              else {
                p_raw->iBSD[p_raw->nBSD] = bsd_adcs->adc[i].channel;
                p_raw->aBSD[p_raw->nBSD] = bsd_adcs->adc[i].value;
                p_raw->tBSD[p_raw->nBSD] = ntdcs;
                p_raw->nBSD++;
              }
            }
          }
        }
  
        p_raw->nBGV = 0;
        if ((bgv_adcs=data_getGroup(event,GROUP_BGV_ADCS,0)) != NULL) {
          for (i=0;i<bgv_adcs->nadc;i++) {
            if (bgv_adcs->adc[i].value >= AMIN_BGV) {
              if (p_raw->nBGV > NT_MAX_BGV_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in bgv, maximum allowed was %d\n",
                          p_raw->nBGV, NT_MAX_BGV_HIT);
                break;
              }
              else {
                p_raw->iBGV[p_raw->nBGV] = bgv_adcs->adc[i].channel;
                p_raw->aBGV[p_raw->nBGV] = bgv_adcs->adc[i].value;
                p_raw->tBGV[p_raw->nBGV] = 0;
                p_raw->nBGV++;
              }
            }
          }
        }
        if ((bgv_tdcs=data_getGroup(event,GROUP_BGV_TDCS,0)) != NULL) {
          for (i=0;i<bgv_tdcs->ntdc;i++) {
            if ((bgv_tdcs->tdc[i].le >= TMIN_BGV) &&
                (bgv_tdcs->tdc[i].le <= TMAX_BGV)) {
              for (j=0;j<p_raw->nBGV;j++) {
                if (p_raw->iBGV[j] == bgv_tdcs->tdc[i].channel)
                  break;
              }
              if (j < p_raw->nBGV)
                p_raw->tBGV[j] = bgv_tdcs->tdc[i].le;
              else if (p_raw->nBGV > NT_MAX_BGV_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in bgv, maximum allowed was %d\n",
                          p_raw->nBGV, NT_MAX_BGV_HIT);
                break;
              }
              else {
                p_raw->iBGV[p_raw->nBGV] = bgv_tdcs->tdc[i].channel;
                p_raw->tBGV[p_raw->nBGV] = bgv_tdcs->tdc[i].le;
                p_raw->aBGV[p_raw->nBGV] = 0;
                p_raw->nBGV++;
              }
            }
          }
        }
  
        p_raw->nCPV = 0;
        if ((cpv_tdcs=data_getGroup(event,GROUP_CPV_TDCS,0)) != NULL) {
          for (i=0;i<cpv_tdcs->ntdc;i++) {
            if ((cpv_tdcs->tdc[i].le >= TMIN_CPV) &&
                (cpv_tdcs->tdc[i].le <= TMAX_CPV)) {
              if (p_raw->nCPV > NT_MAX_CPV_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in cpv, maximum allowed was %d\n",
                          p_raw->nCPV, NT_MAX_CPV_HIT);
                break;
              }
              else {
                p_raw->iCPV[p_raw->nCPV] = cpv_tdcs->tdc[i].channel;
                p_raw->tCPV[p_raw->nCPV] = cpv_tdcs->tdc[i].le;
                p_raw->aCPV[p_raw->nCPV] = 0;
                p_raw->nCPV++;
              }
            }
          }
        }
        if ((cpv_adcs=data_getGroup(event,GROUP_CPV_ADCS,0)) != NULL) {
          for (i=0;i<cpv_adcs->nadc;i++) {
            if (cpv_adcs->adc[i].value >= AMIN_CPV) {
              int ntdcs=0;
              for (j=0;j<p_raw->nCPV;j++) {
                if (p_raw->iCPV[j] == cpv_adcs->adc[i].channel) {
                  p_raw->aCPV[j] = cpv_adcs->adc[i].value;
                  ntdcs++;
                }
              }
              if (p_raw->nCPV > NT_MAX_CPV_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in cpv, maximum allowed was %d\n",
                          p_raw->nCPV, NT_MAX_CPV_HIT);
                break;
              }
              else {
                p_raw->iCPV[p_raw->nCPV] = cpv_adcs->adc[i].channel;
                p_raw->aCPV[p_raw->nCPV] = cpv_adcs->adc[i].value;
                p_raw->tCPV[p_raw->nCPV] = ntdcs;
                p_raw->nCPV++;
              }
            }
          }
        }
  
        p_raw->nUPV = 0;
        if ((upv_tdcs=data_getGroup(event,GROUP_UPV_TDCS,0)) != NULL) {
          for (i=0;i<upv_tdcs->ntdc;i++) {
            if ((upv_tdcs->tdc[i].le >= TMIN_UPV) &&
                (upv_tdcs->tdc[i].le <= TMAX_UPV)) {
              if (p_raw->nUPV > NT_MAX_UPV_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in upv, maximum allowed was %d\n",
                          p_raw->nUPV, NT_MAX_UPV_HIT);
                break;
              }
              else {
                p_raw->iUPV[p_raw->nUPV] = upv_tdcs->tdc[i].channel;
                p_raw->tUPV[p_raw->nUPV] = upv_tdcs->tdc[i].le;
                p_raw->aUPV[p_raw->nUPV] = 0;
                p_raw->nUPV++;
              }
            }
          }
        }
        if ((upv_adcs=data_getGroup(event,GROUP_UPV_ADCS,0)) != NULL) {
          for (i=0;i<upv_adcs->nadc;i++) {
            if (upv_adcs->adc[i].value >= AMIN_UPV) {
              int ntdcs=0;
              for (j=0;j<p_raw->nUPV;j++) {
                if (p_raw->iUPV[j] == upv_adcs->adc[i].channel) {
                  p_raw->aUPV[j] = upv_adcs->adc[i].value;
                  ntdcs++;
                }
              }
              if (p_raw->nUPV > NT_MAX_UPV_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in upv, maximum allowed was %d\n",
                          p_raw->nUPV,NT_MAX_UPV_HIT);
                break;
              }
              else {
                p_raw->iUPV[p_raw->nUPV] = upv_adcs->adc[i].channel;
                p_raw->aUPV[p_raw->nUPV] = upv_adcs->adc[i].value;
                p_raw->tUPV[p_raw->nUPV] = ntdcs;
                p_raw->nUPV++;
              }
            }
          }
        }
  
        p_raw->nTAG = 0;
        if ((tagger_tdcs=data_getGroup(event,GROUP_TAGGER_TDCS,0)) != NULL) {
          for (i=0;i<tagger_tdcs->ntdc;i++) {
            if ((tagger_tdcs->tdc[i].le >= TMIN_TAG) &&
                (tagger_tdcs->tdc[i].le <= TMAX_TAG)) {
              if (p_raw->nTAG > NT_MAX_TAG_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in tagger, maximum allowed was %d\n",
                          p_raw->nTAG,NT_MAX_TAG_HIT);
                break;
              }
              else {
                p_raw->iTAG[p_raw->nTAG] = tagger_tdcs->tdc[i].channel;
                p_raw->tTAG[p_raw->nTAG] = tagger_tdcs->tdc[i].le;
                p_raw->nTAG++;
              }
            }
          }
        }
  
        p_raw->nLGD = 0;
        if ((lgd_adcs=data_getGroup(event,GROUP_LGD_ADCS,0)) != 0) {
          for (i=0;i<lgd_adcs->nadc;i++) {
            if (lgd_adcs->adc[i].value >= AMIN_LGD) {
              if (p_raw->nLGD > NT_MAX_LGD_HIT) {
                fprintf(stderr,"Warning - overflow of ntuple arrays: ");
                fprintf(stderr,"%d hits in lgd, maximum allowed was %d\n",
                          p_raw->nLGD,NT_MAX_LGD_HIT);
                break;
              }
              else {
                p_raw->iLGD[p_raw->nLGD] = lgd_adcs->adc[i].channel % 28;
                p_raw->jLGD[p_raw->nLGD] = lgd_adcs->adc[i].channel / 28;
                p_raw->aLGD[p_raw->nLGD] = lgd_adcs->adc[i].value;
                p_raw->nLGD++;
              }
            }
          }
        }
      }

      if (p_bsd || p_pix || p_rec || p_clust || p_coin || p_talign || p_mes) {
        make_bsd_hits(event,NULL,NT_MAX_BSD_HIT);
      }
      if (p_bgv || p_clust || p_mes) {
        make_bgv_hits(event,NULL,NT_MAX_BGV_HIT);
      }
      if (p_cpv || p_coin) {
        make_cpv_hits(event,NULL,NT_MAX_CPV_HIT);
      }
      if (p_upv) {
        make_upv_hits(event,NULL,NT_MAX_UPV_HIT);
      }
      if (p_tag || p_coin || p_talign) {
        make_tagger_hits(event,NULL,NT_MAX_TAG_HIT);
      }
      if (p_lgd || p_clust || p_mes) {
        make_lgd_hits(event,NULL,NT_MAX_LGD_HIT);
      }

      suppress_dead_channels(event);

      if (p_bsd || p_pix || p_bgv || p_cpv || p_upv || p_talign ||
          p_rec || p_tag || p_clust || p_coin || p_mes) {
        make_time_list(event,NULL,NT_MAX_TIME_LIST);
      }
      if (p_pix || p_clust || p_mes || p_coin || p_rec || p_talign) {
        make_bsd_pixels(event,NULL,NT_MAX_BSD_PIX);
      }
      if (p_clust || p_mes || p_rec) {
        make_bgv_clusters(event,NULL,NT_MAX_BCLUSTERS);
      }
      if (p_rec || p_clust || p_talign) {
        make_recoil(event,NULL,NT_MAX_RECOILS);
      }
      if (p_clust || p_mes) {
        make_lgd_clusters(event,NULL,100);
        make_photons(event,NULL,NT_MAX_PHOTONS);
      }
      if (p_coin) {
        make_tagger_clusters(event,NULL,NT_MAX_TAG_HIT);
        make_tagger_photons(event,NULL,NT_MAX_TAG_COIN);
      }
      if (p_mes) {
        make_mesons(event,NULL,NT_MAX_MESONS);
      }

      if (p_bsd) {
        if ((bsd_hits=data_getGroup(event,GROUP_BSD_HITS,0)) != NULL) {
          p_bsd->nhbsd=bsd_hits->nhits;
          for (i=0;i<p_bsd->nhbsd;i++) {
            p_bsd->chbsd[i]=bsd_hits->hit[i].channel;
            p_bsd->Ebsd[i]=bsd_hits->hit[i].energy;
            p_bsd->ntbsd[i]=bsd_hits->hit[i].times;
            p_bsd->t1bsd[i]=bsd_hits->hit[i].tindex;
          }
        }
        else {
          p_bsd->nhbsd=0;
        }
      }

      if (p_pix) {
        if ((bsd_pixs=data_getGroup(event,GROUP_BSD_PIXELS,0)) != NULL) {
          p_pix->npix=bsd_pixs->npixels;
          for (i=0;i<p_pix->npix;i++) {
            p_pix->lpix[i]=bsd_pixs->pixel[i].left;
            p_pix->rpix[i]=bsd_pixs->pixel[i].right;
            p_pix->spix[i]=bsd_pixs->pixel[i].straight;
            p_pix->ipix[i]=bsd_pixs->pixel[i].ring;
            p_pix->zpix[i]=bsd_pixs->pixel[i].z;
            p_pix->phipix[i]=bsd_pixs->pixel[i].phi;
            p_pix->rpixt[i]=bsd_pixs->pixel[i].time[0];
            p_pix->lpixt[i]=bsd_pixs->pixel[i].time[1];
            p_pix->spixt[i]=bsd_pixs->pixel[i].time[2];
            p_pix->rpixE[i]=bsd_pixs->pixel[i].energy[0];
            p_pix->lpixE[i]=bsd_pixs->pixel[i].energy[1];
            p_pix->spixE[i]=bsd_pixs->pixel[i].energy[2];
          }
        }
        else {
          p_pix->npix=0;
        }
      }

      if (p_rec) {
        if ((recoils=data_getGroup(event,GROUP_RECOIL,0)) != NULL) {
          p_rec->nrec=recoils->nrecoils;
          p_rec->trec0=recoils->best_t0;
          for (i=0;i<p_rec->nrec;i++) {
            p_rec->therec[i]=recoils->recoil[i].theta;
            p_rec->phirec[i]=recoils->recoil[i].phi;
            p_rec->dErec[i]=recoils->recoil[i].dEdx;
            p_rec->Erec[i]=recoils->recoil[i].energy;
            p_rec->trec[i]=recoils->recoil[i].time;
            p_rec->mrec[i]=recoils->recoil[i].npixels;
          }
        }
        else {
          p_rec->nrec=0;
        }
      }

      if (p_bclust) {
        if ((bgv_clusters=data_getGroup(event,GROUP_BGV_CLUSTERS,0)) != NULL) {
          p_bclust->nbcl=bgv_clusters->nClusters;
          for (i=0;i<p_bclust->nbcl;i++) {
            p_bclust->bcE[i]=bgv_clusters->clusters[i].energy;
            p_bclust->bcphi[i]=bgv_clusters->clusters[i].phi;
            p_bclust->bcz[i]=bgv_clusters->clusters[i].z;
            p_bclust->bct[i]=bgv_clusters->clusters[i].t;
            p_bclust->bcsE[i]=bgv_clusters->clusters[i].scintEnergy;
          }
        }
        else {
          p_bclust->nbcl=0;
        }
      }

      if (p_bgv) {
        time_list=data_getGroup(event,GROUP_TIME_LIST,0);
        if ((bgv_hits=data_getGroup(event,GROUP_BGV_HITS,0)) != NULL) {
          p_bgv->nhbgv=bgv_hits->nhits;
          for (i=0;i<p_bgv->nhbgv;i++) {
            p_bgv->chbgv[i]=bgv_hits->hit[i].channel;
            p_bgv->EbgvDwn[i]=bgv_hits->hit[i].energy[1];
            p_bgv->EbgvUp[i]=bgv_hits->hit[i].energy[0];
            p_bgv->ntbgvDwn[i]=bgv_hits->hit[i].times[1];
            p_bgv->t1bgvDwn[i]=bgv_hits->hit[i].tindex[1];
            p_bgv->ntbgvUp[i]=bgv_hits->hit[i].times[0];
            p_bgv->t1bgvUp[i]=bgv_hits->hit[i].tindex[0];
          }
        }
        else {
          p_bgv->nhbgv=0;
        }
      }

      if (p_cpv) {
        if ((cpv_hits=data_getGroup(event,GROUP_CPV_HITS,0)) != NULL) {
          p_cpv->nhcpv=cpv_hits->nhits;
          for (i=0;i<p_cpv->nhcpv;i++) {
            p_cpv->chcpv[i]=cpv_hits->hit[i].channel;
            p_cpv->Ecpv[i]=cpv_hits->hit[i].energy;
            p_cpv->ntcpv[i]=cpv_hits->hit[i].times;
            p_cpv->t1cpv[i]=cpv_hits->hit[i].tindex;
          }
        }
        else {
          p_cpv->nhcpv=0;
        }
      }

      if (p_upv) {
        if ((upv_hits=data_getGroup(event,GROUP_UPV_HITS,0)) != NULL) {
          p_upv->nhupv=upv_hits->nhits;
          for (i=0;i<p_upv->nhupv;i++) {
            p_upv->chupv[i]=upv_hits->hit[i].channel;
            p_upv->Eupv[i]=upv_hits->hit[i].energy;
            p_upv->ntupv[i]=upv_hits->hit[i].times;
            p_upv->t1upv[i]=upv_hits->hit[i].tindex;
          }
        }
        else {
          p_upv->nhupv=0;
        }
      }

      if (p_tag) {
        if ((tagger_hits=data_getGroup(event,GROUP_TAGGER_HITS,0)) != NULL) {
          p_tag->nhtag=tagger_hits->nhits;
          for (i=0;i<p_tag->nhtag;i++) {
            p_tag->chtag[i]=tagger_hits->hit[i].channel;
            p_tag->Etag0[i]=tagger_hits->hit[i].energy[0];
            p_tag->Etag1[i]=tagger_hits->hit[i].energy[1];
            p_tag->nttagL[i]=tagger_hits->hit[i].times[0];
            p_tag->nttagR[i]=tagger_hits->hit[i].times[1];
            p_tag->t1tagL[i]=tagger_hits->hit[i].tindex[0];
            p_tag->t1tagR[i]=tagger_hits->hit[i].tindex[1];
          }
        }
        else {
          p_tag->nhtag=0;
        }
      }

      if (p_tlist) {
        if (time_list=data_getGroup(event,GROUP_TIME_LIST,0)) {
          int ntimes=(time_list->nelements < NT_MAX_TIME_LIST)? 
                      time_list->nelements : NT_MAX_TIME_LIST;
          p_tlist->ntimes=ntimes;
          for (i=0;i<ntimes;i++) {
            float tedge=time_list->element[i].le;
            tedge=(tedge > -999.999)? tedge:-999.999;
            tedge=(tedge < +999.999)? tedge:+999.999;
            p_tlist->le[i]=tedge;
          }
        }
        else {
          p_tlist->ntimes=0;
        }
      }

      if (p_coin) {
        if ((tagging=data_getGroup(event,GROUP_TAGGER_PHOTONS,0)) != NULL) {
          p_coin->ncoin=tagging->nphotons;
          for (i=0;i<p_coin->ncoin;i++) {
            p_coin->cochan[i]=tagging->photon[i].channel;
            p_coin->cotime[i]=tagging->photon[i].time;
            p_coin->coenergy[i]=tagging->photon[i].energy;
            p_coin->tagweight[i]=tagging->photon[i].tagging_weight;
          }
        }
        else {
          p_coin->ncoin=0;
        }
      }

      if (p_talign) {
        if ((recoils=data_getGroup(event,GROUP_RECOIL,0)) != NULL) {
          p_talign->taref=recoils->best_t0;
        }
        else {
          p_talign->taref=0;
        }
        if (((tagger_hits=data_getGroup(event,GROUP_TAGGER_HITS,0)) != NULL) &&
            ((time_list=data_getGroup(event,GROUP_TIME_LIST,0)) != NULL)) {
          int ntal=0;
          int i,j,k;
          for (i=0;i<tagger_hits->nhits;i++) {
            k=tagger_hits->hit[i].tindex[0];
            if (ntal+tagger_hits->hit[i].times[0] < NT_MAX_TAG_ALIGN) {
              for (j=0;j<tagger_hits->hit[i].times[0];j++,k++) {
                if (fabs(time_list->element[k].le) < 80) {
                  p_talign->chtal[ntal]=tagger_hits->hit[i].channel;
                  p_talign->tal[ntal++]=time_list->element[k].le;
                }
              }
            }
            k=tagger_hits->hit[i].tindex[1];
            if (ntal+tagger_hits->hit[i].times[1] < NT_MAX_TAG_ALIGN) {
              for (j=0;j<tagger_hits->hit[i].times[1];j++,k++) {
                if (fabs(time_list->element[k].le) < 80) {
                  p_talign->chtal[ntal]=tagger_hits->hit[i].channel+19;
                  p_talign->tal[ntal++]=time_list->element[k].le;
                }
              }
            }
          }
          p_talign->ntal=ntal;
        }
      }

      if (p_lgd) {
        if ((lgd_hits=data_getGroup(event,GROUP_LGD_HITS,0)) != NULL) {
          p_lgd->nhlgd=lgd_hits->nhits;
          for (i=0;i<p_lgd->nhlgd;i++) {
            p_lgd->chlgd[i]=lgd_hits->hit[i].channel;
            p_lgd->Elgd[i]=lgd_hits->hit[i].energy;
            p_lgd->clust[i]=0;
          }
          if ((lgd_hits=data_getGroup(event,GROUP_LGD_CLUSTER_HITS,0)) &&
              (lgd_clusters=data_getGroup(event,GROUP_LGD_CLUSTERS,0))) {
            int ic;
            for (ic=0;ic<lgd_clusters->nClusters;ic++) {
              int i1=lgd_clusters->cluster[ic].firstBlock;
              int i2=lgd_clusters->cluster[ic].nBlocks+i1;
              int ii;
              for (ii=i1;ii<i2;ii++) {
                for (i=0;i<p_lgd->nhlgd;i++) {
                  if (p_lgd->chlgd[i]==lgd_hits->hit[ii].channel) {
                    p_lgd->clust[i]=ic+1;
                  }
                }
              }
            }
          }
        }
        else {
          p_lgd->nhlgd=0;
        }
      }

      if (p_clust) {
        if ((photons=data_getGroup(event,GROUP_PHOTONS,0)) != NULL) {
          while (photons->nPhotons > 3) {
              int count = photons->nPhotons;
              merge_lowmass_cluster_pairs(photons);
              if (photons->nPhotons == count || photons->nPhotons == 3)
                 break;
          }
          p_clust->nphot=photons->nPhotons;
          p_clust->nfrwd=0;
          for (i=0;i<p_clust->nphot;i++) {
            p_clust->pvect[i].t=photons->photon[i].energy;
            p_clust->pvect[i].space.x=photons->photon[i].momentum.x;
            p_clust->pvect[i].space.y=photons->photon[i].momentum.y;
            p_clust->pvect[i].space.z=photons->photon[i].momentum.z;
            if (photons->photon[i].cluster >= 0) {
              ++(p_clust->nfrwd);
            }
          }
        }
        else {
          p_clust->nphot=0;
          p_clust->nfrwd=0;
        }
        p_clust->Efrwd=0;
        if ((lgd_hits=data_getGroup(event,GROUP_LGD_HITS,0)) != NULL) {
          for (i=0;i<lgd_hits->nhits;i++) {
            p_clust->Efrwd+=lgd_hits->hit[i].energy;
          }
        }
      }

      if (p_mes) {
        if ((mesons=data_getGroup(event,GROUP_MESONS,0)) != NULL) {
          p_mes->nmes=0;
          for (i=0;i<mesons->nMesons;i++) {
            p_mes->mtype[p_mes->nmes]=mesons->meson[i].type; 
            p_mes->ptot[p_mes->nmes].t=mesons->meson[i].energy;
            p_mes->ptot[p_mes->nmes].space.x=mesons->meson[i].momentum.x;
            p_mes->ptot[p_mes->nmes].space.y=mesons->meson[i].momentum.y;
            p_mes->ptot[p_mes->nmes].space.z=mesons->meson[i].momentum.z;
            p_mes->amass[p_mes->nmes]=mesons->meson[i].mass; 
            p_mes->idtype[p_mes->nmes][0]=mesons->meson[i].dtype[0]; 
            p_mes->idtype[p_mes->nmes][1]=mesons->meson[i].dtype[1]; 
            p_mes->ichild[p_mes->nmes][0]=mesons->meson[i].daughter[0]; 
            p_mes->ichild[p_mes->nmes][1]=mesons->meson[i].daughter[1]; 
            p_mes->nmes++;
          }
        }
        else {
          p_mes->nmes=0;
        }
      }

#if defined HISTOGRAM_ONLY      

/* For run when just histograms are required */

# if defined TAGGER_CALIBRATION
      if (recoils->best_t0 != 0 && recoils->nrecoils == 1) {
        float t0=recoils->best_t0;
        int ice;
        if (hexist_(ntid) == 0) {
          for (ice=0; ice < 38; ice++) {
            hbook2(*ntid+ice,"tagger time",100,-25.,25.,1400,7400.,8800.,0);
          }
        }
        for (i=0; i < tagger_hits->nhits; i++) {
          int il,ir,it;
          ice=tagger_hits->hit[i].channel;
          it=tagger_hits->hit[i].tindex[0];
          for (il=0; il < tagger_hits->hit[i].times[0]; il++,it++) {
            float t=time_list->element[it].le;
            hfill(*ntid+ice,t-t0,(float)event->runNo,1.);
          }
          ice+=19;
          it=tagger_hits->hit[i].tindex[1];
          for (ir=0; ir < tagger_hits->hit[i].times[1]; ir++,it++) {
            float t=time_list->element[it].le;
            hfill(*ntid+ice,t-t0,(float)event->runNo,1.);
          }
        }
      }
# else
      if (hexist_(ntid) == 0) {
        int cntr, ring;
        for (cntr=0; cntr < 48; cntr++) {
          for (ring=0; ring < 1; ring++) {
            hbook2(*ntid+ring*100+cntr,"tagger-pixel time",100,-25.,25.,1400,7400.,8800.,0);
          }
        }
      }
      if (bsd_pixs->npixels == 1) {
        for (i=0; i < tagger_hits->nhits; i++) {
          if (tagger_hits->hit[i].channel == 0) {
            int it=tagger_hits->hit[i].tindex[0];
            int il;
            for (il=0; il < tagger_hits->hit[i].times[0]; il++,it++) {
              int cntr, ring=0;//bsd_pixs->pixel[0].ring;
              float ttag, tpix;
              ttag=time_list->element[it].le;
              cntr=bsd_pixs->pixel[0].right;
              tpix=bsd_pixs->pixel[0].time[0];
              hfill(*ntid+ring*100+cntr,ttag-tpix,(float)event->runNo,1.);
              cntr=bsd_pixs->pixel[0].left;
              tpix=bsd_pixs->pixel[0].time[1];
              hfill(*ntid+ring*100+cntr,ttag-tpix,(float)event->runNo,1.);
              cntr=bsd_pixs->pixel[0].straight;
              tpix=bsd_pixs->pixel[0].time[2];
              hfill(*ntid+ring*100+cntr,ttag-tpix,(float)event->runNo,1.);
            }
          }
        }
      }
# endif
#else

      cleanlgd = lgd_cluster_cleanup(event);
      if (p_raw)
        p_raw->evType += (cleanlgd << 8);
  #if defined CLUSTER_CLEANUP                      
      if (cleanlgd == 0)
  #endif                      
      {
        if (hexist_(ntid) == 0) {
          declare_raw_ntuple(*ntid,p_head);
          switch (event->eventType) {
            case EV_DATA:
              if (p_raw) appendnt_raw_event(*ntid,p_raw);
              if (p_bsd) appendnt_bsd_hits(*ntid,p_bsd);
              if (p_pix) appendnt_bsd_pixs(*ntid,p_pix);
              if (p_bgv) appendnt_bgv_hits(*ntid,p_bgv);
              if (p_cpv) appendnt_cpv_hits(*ntid,p_cpv);
              if (p_upv) appendnt_upv_hits(*ntid,p_upv);
              if (p_tag) appendnt_tag_hits(*ntid,p_tag);
              if (p_tlist) appendnt_time_list(*ntid,p_tlist);
              if (p_rec) appendnt_recoil(*ntid,p_rec);
              if (p_coin) appendnt_tag_coin(*ntid,p_coin);
              if (p_talign) appendnt_tag_align(*ntid,p_talign);
              if (p_lgd) appendnt_lgd_hits(*ntid,p_lgd);
              if (p_clust) appendnt_photons(*ntid,p_clust);
              if (p_bclust) appendnt_bclusters(*ntid,p_bclust);
              if (p_mes) appendnt_mesons(*ntid,p_mes);
            case EV_MCGEN:
              break;
          }
        }
#if defined FULLY_CONTAINED_FORWARD
        {
          int neutral_bgv_hits=0;
          for (i=0; i<p_bgv->nhbgv; i++) {
            ++neutral_bgv_hits;
            for (j=0; j<p_pix->npix; j++) {
              int s1=(p_pix->spix[j]-22)%24;
              int s2=(p_pix->spix[j]-23)%24;
              if (p_bgv->chbgv[i] == s1 ||
                  p_bgv->chbgv[i] == s2) {
                --neutral_bgv_hits;
                break;
              }
            }
          }
          if (p_clust->nphot >= 1 && p_clust->nphot <= 7
              && neutral_bgv_hits == 0
              && p_bgv->nhbgv <= 1
              && p_rec->nrec == 1
              )
#  define FILL_NTUPLE 0
#  if FILL_NTUPLE
            hfnt_(ntid);
#  else
          {
             int id=100;
             int m;
             float ptot[4];
             float masstot;
             float showmult;
             if (hexist_(&id) == 0) {
                char title[80];
                for (id=0; id < 19; id++) {
                   sprintf(title,"tag time channel %d",id);
                   hbook1(id+100,title,200,-50.,50.,0);
                   sprintf(title,"tag cpvbar coinc channel %d",id);
                   hbook1(id+1100,title,200,-50.,50.,0);
                   sprintf(title,"tag cpvbar accid channel %d",id);
                   hbook1(id+2100,title,200,-50.,50.,0);
                   for (m=2; m < 8; m++) {
                      sprintf(title,"tag time channel %d",id);
                      hbook1(id+m*100,title,200,-50.,50.,0);
                      sprintf(title,"tag cpvbar coinc channel %d",id);
                      hbook1(id+m*100+1000,title,200,-50.,50.,0);
                      sprintf(title,"tag cpvbar accid channel %d",id);
                      hbook1(id+m*100+2000,title,200,-50.,50.,0);
                      sprintf(title,"energy channel %d",id);
                      hbook1(id+m*100+20,title,200,0.,10.,0);
                      sprintf(title,"energy coinc channel %d",id);
                      hbook1(id+m*100+1020,title,200,0.,10.,0);
                      sprintf(title,"energy accid channel %d",id);
                      hbook1(id+m*100+2020,title,200,0.,10.,0);
                      sprintf(title,"mass channel %d",id);
                      hbook1(id+m*100+40,title,200,0.,2.,0);
                      sprintf(title,"mass coinc channel %d",id);
                      hbook1(id+m*100+1040,title,200,0.,2.,0);
                      sprintf(title,"mass accid channel %d",id);
                      hbook1(id+m*100+2040,title,200,0.,2.,0);
                   }
                }
                sprintf(title,"LGD shower multiplicity");
                hbook1(90,title,10,0.,10.,0);
             }
             ptot[0]=ptot[1]=ptot[2]=ptot[3]=0;
             for (i=0; i < p_clust->nphot; i++) {
                ptot[0] += p_clust->pvect[i].t;
                ptot[1] += p_clust->pvect[i].space.x;
                ptot[2] += p_clust->pvect[i].space.y;
                ptot[3] += p_clust->pvect[i].space.z;
             }
             masstot = sqrt(sqr(ptot[0])-
                            sqr(ptot[1])-
                            sqr(ptot[2])-
                            sqr(ptot[3]));
             for (i=0; i < p_coin->ncoin; i++) {
                float dt = p_rec->trec0 - p_coin->cotime[i];
                id = p_coin->cochan[i];
                m = p_clust->nphot;
                hfill(id+100,dt,0.,1.);
                hfill(id+m*100,dt,0.,1.);
                hfill(id+m*100+20,ptot[0],0.,1.);
                hfill(id+m*100+40,masstot,0.,1.);
                if (p_coin->vetot[i] == 0) {
                   hfill(id+1100,dt,0.,1.);
                   hfill(id+m*100+1000,dt,0.,1.);
                   if (p_coin->coint[i] > 0) {
                      hfill(id+m*100+1020,ptot[0],0.,1.);
                      hfill(id+m*100+1040,masstot,0.,1.);
                   }
                }
                if (p_coin->vetoa[i] == 0) {
                   hfill(id+2100,dt,0.,1.);
                   hfill(id+m*100+2000,dt,0.,1.);
                   if (p_coin->coina[i] > 0) {
                      hfill(id+m*100+2020,ptot[0],0.,1.);
                      hfill(id+m*100+2040,masstot,0.,1.);
                   }
                }
             }
             showmult=p_clust->nphot;
             hfill(90,showmult,0.,1.);
          }
#  endif
        }
#else
        hfnt_(ntid);
#endif
      }
#endif
      ++srun[1];
      
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
      fprintf(stderr,"Event type : %d\n",event->eventType);
    }          

    if (quest_.iq[0]) {
      ret = OUTPUT_FILE_OVERRUN;
      break;
    }

    for (i=10000; i < records; i*=10) {}
    i /= 10;
    if (records/i*i == records)
      printf("%d\n",records);
    if (records/10000*10000 == records) {
      char chpath[80];
      int nlevel=99;
      hcdir_(chpath,"R",80,1);
      rzstat_(chpath,&nlevel,"Q",80,1);
      if (quest_.iq[11] > 480000000) {
        ret = OUTPUT_FILE_FULL;
        break;
      }
    }
  }

  switch (ret) {
  case DATAIO_BADCRC:
  case DATAIO_ERROR:
    fprintf(stderr,"Error occurred in event %d after %d bytes read.\n",
            records,bytes);
  case OUTPUT_FILE_FULL:
    break;
  default:
    evt_close(fd);
    fd=0;
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
  fprintf(stderr,"%d total bytes read\n",bytes);

  return ret;
}
