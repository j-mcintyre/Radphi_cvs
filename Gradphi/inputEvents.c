/* fetches MonteCarlo event stubs from an event generator file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <disData.h>
#include <dataIO.h>

#include <mctypes.h>

#define BEAM_TARGET_X_OFFSET_CM 0.5
#define BEAM_TARGET_Y_OFFSET_CM 0.0

#define BUFSIZE 100000   /* Something really big */

void gsvertx_(float vert[3], int *ntbeam, int *nttarg, float *tofg,
              float *ubuf, int *nubuf, int *nvtx);
void gskine_(float plab[3], int *ipart, int *nv, float *ubuf, 
             int *nubuf, int *nt);
void gsveru_(int *nvtx, int *nubuf, float *ubuf, int *iadr);
void gskinu_(int *nt, int *nubuf, float *ubuf, int *iadr);
void gsrung_(int *nubuf, float *ubuf, int* iadr);

int loadKinematics(itape_header_t *event)
{
  float *vert;
  int ntbeam;
  int nttarg;
  float ubuf[1];
  int nubuf;
  int nvtx;
  float *plab;
  float *tofg;
  int ipart;
  int nt;
  int iv, ip;
  int numed;
  
  mc_event_t *mc_event;
  mc_vertex_t *mc_vertex;
  mc_particle_t *particle;

  if ((mc_event=data_getGroup(event,GROUP_MC_EVENT,0)) == NULL) {
    fprintf(stderr,"ERROR - MC event group missing from input data\n");
    return 0;
  }

  if (mc_event->nfinalVert > 0) {
    for (iv=1;iv<=mc_event->nfinalVert;iv++) {
      if ((mc_vertex=data_getGroup(event,GROUP_MC_EVENT,iv)) == NULL) {
        fprintf(stderr,"ERROR - MC event subgroup missing from input data\n");
        return 0;
      }
      /* add an offset to simulate misalignment of beam during running */
      mc_vertex->origin.space.x += BEAM_TARGET_X_OFFSET_CM;
      mc_vertex->origin.space.y += BEAM_TARGET_Y_OFFSET_CM;
      vert=&mc_vertex->origin.space.x;
      tofg=&mc_vertex->origin.t;

/* require vertices to lie inside the target in order to be tracked */
#define MED_BERYLLIUM_TGT 22    /* must match value in geometry.inc */
      gmedia_(vert,&numed);

      ntbeam=0;
      nttarg=0;
      nubuf=sizeof_mc_vertex_t(mc_vertex->npart)/4;
      gsvertx_(vert,&ntbeam,&nttarg,tofg,(float *)mc_vertex,&nubuf,&nvtx);
      if (nvtx == 0) {
        fprintf(stderr,"WARNING - gsvert refused to store vertex!\n");
        break;
      }
      else if (numed != MED_BERYLLIUM_TGT) {
        break;
      }
      else if (mc_vertex->npart > 0) {
        for (ip=0;ip<mc_vertex->npart;ip++) {
          particle=&mc_vertex->part[ip];
          plab=&particle->momentum.space.x;
          ipart=particle->kind;
          nubuf=0;
          gskine_(plab,&ipart,&nvtx,ubuf,&nubuf,&nt);
          if (nt == 0) {
            fprintf(stderr,"WARNING - gskine refused to store track!\n");
            break;
          }
        }
      }
    }
  }

/* attach the mc_event group to the run structure for output later */

  {
    int iadr=0;
    nubuf=sizeof_mc_event_t(mc_event->npart)/4;
    gsrung_(&nubuf,(float *)mc_event,&iadr);
  }
  return 1;
}

int store_mc_event_(mc_event_t *ubuf)
{
   ubuf->npart = 0;
   ubuf->type = 1;
   ubuf->weight = 1;
   ubuf->nfinalVert = 1;
   return sizeof_mc_event_t(0);
}

int store_mc_vertex_(mc_vertex_t *ubuf, int *kind, float vtx[3],
	                        	float mom[3], float *Etot)
{
   ubuf->npart = 1;
   ubuf->kind = 0;
   ubuf->origin.space.x = vtx[0];
   ubuf->origin.space.y = vtx[1];
   ubuf->origin.space.z = vtx[2];
   ubuf->part[0].kind = *kind;
   ubuf->part[0].momentum.space.x = mom[0];
   ubuf->part[0].momentum.space.y = mom[1];
   ubuf->part[0].momentum.space.z = mom[2];
   ubuf->part[0].momentum.t = *Etot;
   return sizeof_mc_vertex_t(1);
}

int load_next_event_(char *infile, int len)
{
  static FILE *fp=NULL;
  static itape_header_t *event=NULL;
  static int howMany, firstOne;
  char filename[120];
  int specs;
  int ret;
  int i;
  

  /* Malloc the event buffer - only do this once */
  if ((event == NULL) && ((event=malloc(BUFSIZE)) == NULL)) {
    fprintf(stderr,"ERROR - failed to allocating input event buffer!\n");
    fclose(fp);
    fp=NULL;
    return 0;
  }
  
  if (fp == NULL) {
    specs=sscanf(infile,"%s %d %d",&filename,&firstOne,&howMany);
    if (specs < 3) {
      howMany=999999999;
    }
    if (specs < 2) {
      firstOne=1;
    }
    if (specs < 1) {
      fprintf(stderr,"WARNING - invalid input event filename\n");
      return 0;
    }
    if ((fp=fopen(filename,"r")) == NULL) {
      fprintf(stderr,"WARNING - failed to open input event file %s\n",filename);
      return 0;
    }
    else {
      printf("INFO - starting simulation on input file %s\n",filename);
    }
    if (firstOne > 1) {
      printf("INFO - skipping to start at event %d\n",firstOne);
      while (--firstOne > 0) {
        if ((ret=data_read(fileno(fp),event,BUFSIZE)) != DATAIO_OK) {
          fclose(fp);
          fp=NULL;
          return 0;
        }
      }
    }
  }
  
  /* Note that data_read returns DATAIO_EOT at the EOD on tapes
     and at the EOF on files */

  if (howMany-- == 0 ||
     ((ret=data_read(fileno(fp),event,BUFSIZE)) != DATAIO_OK)) {
    fclose(fp);
    fp=NULL;
    return 0;
  }

  return loadKinematics(event);
}
