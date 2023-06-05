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
#include <cernlib.h>
#include <mctypes.h>

#include "mctuple.h"
#include "hituple.h"

//#define HISTOGRAM_ONLY 1
//#define CLUSTER_CLEANUP 1
//#define FULLY_CONTAINED_FORWARD 1

#define sqr(x) ((x)*(x))
#define OUTPUT_FILE_FULL 95
#define OUTPUT_FILE_OVERRUN 96

void merge_lowmass_cluster_pairs(photons_t *photons);

extern struct {
  int iq[100];
} quest_;

int fd=0;

int mctuple_ (char *inspec, int *ntid, int *count, int *start, int len)
{
  int ret=0;
  int srun[]={0,0};
  static itape_header_t *event=NULL;
  int runNo=999999;
  char infile[1000];
  char *fin;
  int i,j;
  
  static nt_head_t head;
  static nt_head_t *p_head=&head;
  static nt_mc_event_t *p_mc=NULL;
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

  mc_event_t *mc_event=NULL;
  mc_vertex_t *mc_vertex=NULL;
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

  int bytes=0;
  int records=0;

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
    if (p_mc)    free(p_mc), p_mc = NULL;
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
      sprintf(tok,"garbxvcutlpmoe");
    }
    if (index(tok,'g') != NULL) {
      p_mc = malloc(sizeof(nt_mc_event_t));
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
    case EV_MCGEN:
      if (event->runNo != srun[0]) {
        if (setup_makeHits(event)) {
          fprintf(stderr,
                  "Error in mctuple: setup_makeHits failed, quitting!\n");
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

      if (p_mc) {
        if (mc_event=data_getGroup(event,GROUP_MC_EVENT,0)) {
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
        }
        else {
          p_mc->nI=0;
          p_mc->type=0;
          p_mc->wgt=0.;
          p_mc->nV=0;
          p_mc->nF=0;
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
      if (hexist_(ntid) == 0) {
        hbook1(*ntid,"photon multiplicity",20,-0.5,19.5,0);
      }     
      if (photons) {
        hfill(*ntid,(float)photons->nPhotons,0.,1.);
      }
#else

      if (p_clust // && (p_clust->nphot >= 2) 
         )
      {
        if (hexist_(ntid) == 0) {
          declare_mc_ntuple(*ntid,p_head);
          switch (event->eventType) {
            case EV_DATA:
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
              if (p_mc) appendnt_mc_event(*ntid,p_mc);
          }
        }

        p_mc->type = lgd_cluster_cleanup(event);
# if defined CLUSTER_CLEANUP   
        if (p_mc->type == 0)
# endif

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
          if (neutral_bgv_hits == 0 &&
              p_bgv->nhbgv <= 1 &&
              p_rec->nrec == 1 &&
              p_clust->nphot >= 1 &&
              p_clust->nphot <= 7) {
            hfnt_(ntid);
          }
        }
#else
        hfnt_(ntid);
#endif
      }
#endif
      ++srun[1];
      break;
    default:
      fprintf(stderr,"Error - incorrect event type: %d !\n",event->eventType);
      fprintf(stderr,"in the event %d\n",records);
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

  fprintf(stderr,"Processed %d records\n",records);
  return ret;
}
