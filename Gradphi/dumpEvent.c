/* fetches MonteCarlo event stubs from an event generator file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <disData.h>
#include <dataIO.h>

#include <mctypes.h>
#include <lgdCluster.h>
#include <iitypes.h>

#define BUFSIZE 100000   /* Something really big */

void Usage(void)
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"  dumpEvent <options> file\n");
  fprintf(stderr,"Options are:\n");
  fprintf(stderr,"\t-h\t\tPrint this message\n");
  fprintf(stderr,"\t-M#\t\tProcess at most # events\n");
}

/* prints the current mc_event group */
int mcEventPrint (mc_event_t *mc_event) {

 int k=0;
 
 printf("Printing group MC_EVENT\n");
 printf(" type: %i\n",mc_event->type);
 printf(" weight: %f\n",mc_event->weight);
 printf(" Printing Initial State\n");
 printf("  Number of particles: %i\n",mc_event->npart);
 for (k=0;k<mc_event->npart;k++) {
  printf("   Particle number %i\n",(k+1));
  printf("   Kind: %i\n",mc_event->part[k].kind);
  printf("   px: %f\n",mc_event->part[k].momentum.space.x);
  printf("   py: %f\n",mc_event->part[k].momentum.space.y);
  printf("   pz: %f\n",mc_event->part[k].momentum.space.z);
  printf(" Printing Final State\n");
  printf("  Number of vertices: %i\n",mc_event->nfinalVert);
 }
 return 1; 
}

/* prints the current mc_vertex group */
int mcVertexPrint (int iv, mc_vertex_t *mc_vertex) {

 int l=0;

 printf("   Vertex number %i\n",iv);
 printf("   Kind: %i\n",mc_vertex->kind);
 printf("   x of Origin: %f\n",mc_vertex->origin.space.x);
 printf("   y of Origin: %f\n",mc_vertex->origin.space.y);
 printf("   z of Origin: %f\n",mc_vertex->origin.space.z);
 printf("   Number of particles in the vertex: %i\n",mc_vertex->npart);
 for (l=0;l<mc_vertex->npart;l++) {
   printf("    Particle number %i\n",(l+1));
   printf("     Kind: %i\n",mc_vertex->part[l].kind);
   printf("     px : %f\n",mc_vertex->part[l].momentum.space.x);
   printf("     py : %f\n",mc_vertex->part[l].momentum.space.y);
   printf("     pz : %f\n",mc_vertex->part[l].momentum.space.z);
 } 
 printf("------------------End of MC_GROUP--------------------------\n");
 return 1; 
}

/* lgd hits */
int lgdHitsPrint(lgd_hits_t *lgd_hits)
 {

 int k=0;
 
 printf("Printing group LGD_HITS\n");
 printf(" Number of hits: %i\n",lgd_hits->nhits);
 for (k=0;k<lgd_hits->nhits;k++) {
  printf("  Hit number %i\n",(k+1));
  printf("   channel: %i\n",lgd_hits->hit[k].channel);
  printf("   energy: %f\n",lgd_hits->hit[k].energy);
 }
  printf("------------------End of LGD_HITS group----------------------\n");
  return 1;
 } 

/* cpv hits */
int cpvHitsPrint(cpv_hits_t *cpv_hits)
 {

 int k=0;
 
 printf("Printing group CPV_HITS\n");
 printf(" Number of hits: %i\n",cpv_hits->nhits);
 for (k=0;k<cpv_hits->nhits;k++) {
  printf("  Hit number %i\n",(k+1));
  printf("   channel: %i\n",cpv_hits->hit[k].channel);
  printf("   energy: %f\n",cpv_hits->hit[k].energy);
  printf("   Number of timehits: %i\n",cpv_hits->hit[k].times);
  printf("   Position in the timelist: %i\n",cpv_hits->hit[k].tindex);
 }
  printf("------------------End of CPV_HITS group----------------------\n");
  return 1;
 } 
 
/* rpd hits */
int rpdHitsPrint(rpd_hits_t *rpd_hits)
 {

 int k=0;
 
 printf("Printing group RPD_HITS\n");
 printf(" Number of hits: %i\n",rpd_hits->nhits);
 for (k=0;k<rpd_hits->nhits;k++) {
  printf("  Hit number %i\n",(k+1));
  printf("   channel: %i\n",rpd_hits->hit[k].channel);
  printf("   energy: %f\n",rpd_hits->hit[k].energy);
  printf("   Number of timehits: %i\n",rpd_hits->hit[k].times);
  printf("   Position in the timelist: %i\n",rpd_hits->hit[k].tindex);
 }
  printf("------------------End of RPD_HITS group----------------------\n");
  return 1;
 }
  
/* upv hits */
int upvHitsPrint(upv_hits_t *upv_hits)
 {

 int k=0;
 
 printf("Printing group UPV_HITS\n");
 printf(" Number of hits: %i\n",upv_hits->nhits);
 for (k=0;k<upv_hits->nhits;k++) {
  printf("  Hit number %i\n",(k+1));
  printf("   channel: %i\n",upv_hits->hit[k].channel);
  printf("   energy: %f\n",upv_hits->hit[k].energy);
  printf("   Number of timehits: %i\n",upv_hits->hit[k].times);
  printf("   Position in the timelist: %i\n",upv_hits->hit[k].tindex);
 }
  printf("------------------End of UPV_HITS group----------------------\n");
  return 1;
 }
  
/* BSD hits */
int bsdHitsPrint(bsd_hits_t *bsd_hits)
 {
 int k=0;
 printf("Printing group BSD_HITS\n");
 printf(" Number of hits: %i\n",bsd_hits->nhits);
 for(k=0;k<bsd_hits->nhits;k++) {
  printf("  Hit number %i\n",(k+1));
  printf("   channel: %i\n",bsd_hits->hit[k].channel);
  printf("   energy: %f\n",bsd_hits->hit[k].energy);
  printf("   Number of timehits: %i\n",bsd_hits->hit[k].times);
  printf("   Position in the timelist: %i\n",bsd_hits->hit[k].tindex);
  }
 printf("------------------End of BSD_HITS group----------------------\n");
 return 1;
 }

/* BGV hits */
int bgvHitsPrint(bgv_hits_t *bgv_hits)
 {
 int k=0;
 printf("Printing group BGV_HITS\n");
 printf(" Number of hits: %i\n",bgv_hits->nhits);
 for(k=0;k<bgv_hits->nhits;k++) {
  printf("  Hit number %i\n",(k+1));
  printf("   channel: %i\n",bgv_hits->hit[k].channel);
  printf("   upstream energy: %f\n",bgv_hits->hit[k].energy[0]);
  printf("   downstream energy: %f\n",bgv_hits->hit[k].energy[1]);
  printf("   Number of upstream timehits: %i\n",bgv_hits->hit[k].times[0]);
  printf("   Position in the timelist: %i\n",bgv_hits->hit[k].tindex[0]);
  printf("   Number of downstream timehits: %i\n",bgv_hits->hit[k].times[1]);
  printf("   Position in the timelist: %i\n",bgv_hits->hit[k].tindex[1]);
  }
 printf("------------------End of BGV_HITS group----------------------\n");
 return 1;
 }

/* time list */
int timeListPrint(time_list_t *time_list)
 {

 int k=0;
 
 printf("Printing group TIME_LIST\n");
 printf(" Number of elements: %i\n",time_list->nelements);
 for (k=0;k<time_list->nelements;k++) {
  printf("  Element number %i\n",(k+1));
  printf("  Time of the hit: %e\n",time_list->element[k].le);
 }
  printf("------------------End of TIME_LIST group----------------------\n");
  return 1;
 }

/* prints the current event */
int printEvent (int noMC,
                mc_event_t * mc_event,
                mc_vertex_t * mc_vertex[],
                lgd_hits_t * lgd_hits,
                cpv_hits_t * cpv_hits,
                rpd_hits_t * rpd_hits,
                upv_hits_t * upv_hits,
                bsd_hits_t * bsd_hits,
		bgv_hits_t * bgv_hits,
                time_list_t *time_list) {

 if (noMC > 0) {
  int iv;
  if (!mcEventPrint(mc_event)) return 0;
  for (iv=0;iv<noMC;iv++) {
    if (!mcVertexPrint(iv+1,mc_vertex[iv])) return 0;
  }
 }
 if (lgd_hits->nhits) if (!lgdHitsPrint(lgd_hits)) return 0;
 if (cpv_hits->nhits) if (!cpvHitsPrint(cpv_hits)) return 0;
/* if (!rpdHitsPrint(rpd_hits)) return 0;*/
 if (upv_hits->nhits) if (!upvHitsPrint(upv_hits)) return 0;
 if (bsd_hits->nhits) if (!bsdHitsPrint(bsd_hits)) return 0;
 if (bgv_hits->nhits) if (!bgvHitsPrint(bgv_hits)) return 0;
 if (time_list->nelements) if (!timeListPrint(time_list)) return 0; 
 return 0;              
}                

/* Gets the next event from the file */
int getEvent(char *filename)
{
  static FILE *fp=NULL;
  static itape_header_t *event=NULL;
  mc_event_t *mc_event=NULL;
  mc_vertex_t *mc_vertex[10];
  lgd_hits_t *lgd_hits=NULL;
  cpv_hits_t *cpv_hits=NULL;
  rpd_hits_t *rpd_hits=NULL;
  upv_hits_t *upv_hits=NULL;
  bsd_hits_t *bsd_hits=NULL;
  bgv_hits_t *bgv_hits=NULL;
  time_list_t *time_list=NULL;
  int ret;
  int i,iv;
  int noMC=0;
  
  if ((fp == NULL) && ((fp=fopen(filename,"r")) == NULL)) {
    fprintf(stderr,"WARNING - failed to open input events file %s\n",filename);
    return 1;
  }
  
  /* Malloc the event buffer - only do this once */
  if ((event == NULL) && ((event=malloc(BUFSIZE)) == NULL)) {
    fprintf(stderr,"ERROR - failed to allocating input event buffer!\n");
    fclose(fp);
    return 1;
  }

  /* Note that data_read returns DATAIO_EOT at the EOD on tapes
     and at the EOF on files */

  if ((ret=data_read(fileno(fp),event,BUFSIZE)) == DATAIO_OK) {
    printf("===>New event, run %d, event %d\n",event->runNo,event->eventNo);
    
    if ((mc_event=data_getGroup(event,GROUP_MC_EVENT,0)) == NULL) {
      fprintf(stderr,"Warning: mc_event group missing event %d\n",event->eventNo);
    }
    else {
      noMC=mc_event->nfinalVert;
    }
    
    for (iv=0;iv<noMC;iv++) {
      if ((mc_vertex[iv]=data_getGroup(event,GROUP_MC_EVENT,iv+1)) == NULL) {
        fprintf(stderr,"Warning: mc_vertex subgroup missing event %d\n",event->eventNo);
      }
    }
    
    if ((lgd_hits=data_getGroup(event,GROUP_LGD_HITS,0)) == NULL) {
     fprintf(stderr,"ERROR - lgd_hits group missing from input data\n");
     fclose(fp);
     return 1;     
    }  
 
    else if ((cpv_hits=data_getGroup(event,GROUP_CPV_HITS,0)) == NULL) {
     fprintf(stderr,"ERROR - cpv_hits group missing from input data\n");
     fclose(fp);
     return 1;
    } 
    else if ((rpd_hits=data_getGroup(event,GROUP_RPD_HITS,0)) == NULL) {
     fprintf(stderr,"ERROR - rpd_hits group missing from input data\n");
     fclose(fp);
     return 1;      
    }
    else if ((upv_hits=data_getGroup(event,GROUP_UPV_HITS,0)) == NULL) {
     fprintf(stderr,"ERROR - upv_hits group missing from input data\n");
     fclose(fp);
     return 1;          
    } 
  
    else if ((bsd_hits=data_getGroup(event,GROUP_BSD_HITS,0)) == NULL) {
     fprintf(stderr,"ERROR - bsd_hits group missing from input data\n");
     fclose(fp);
     return 1;
    }

    else if ((bgv_hits=data_getGroup(event,GROUP_BGV_HITS,0)) == NULL) {
     fprintf(stderr,"ERROR - bgv_group missing from input data\n");
     fclose(fp);
     return 1;
    }

    else if ((time_list=data_getGroup(event,GROUP_TIME_LIST,0)) == NULL) {
     fprintf(stderr,"ERROR - time list group missing from input data\n");
     fclose(fp);
     return 1;          
    } 
    
  }
  else if (ret!=DATAIO_EOT) {
    fprintf(stderr,"ERROR - Couldn't read event group from file %s\n",filename);
    fclose(fp);
    return 1;
  }
  else {
   fclose(fp);
   return 2;
  }

  return printEvent(noMC,mc_event,mc_vertex,lgd_hits,cpv_hits,rpd_hits,upv_hits,
                    bsd_hits,bgv_hits,time_list);
}

/* main program */
int main(int argc, char *argv[])
{
  int i=1;
  char *filename;
  int maxEvents=0;
  int c;
  
  opterr=0;

  while ((c = getopt(argc,argv,"M:h")) != -1)
    switch(c)
      {
      case 'M':
	fprintf(stderr,"maxEv= %s\n",optarg);
	if(optarg!=NULL) maxEvents=atoi(optarg);
	break;
      case 'h':
	Usage();
	exit(EXIT_FAILURE);
      case '?':
	if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        exit(EXIT_FAILURE);
      default:
	Usage();
	exit(EXIT_FAILURE);
      }
  
  if(*argv[optind] != '-'){
    filename=argv[optind];
      while (!getEvent(filename)) {
	printf("==============================================================\n");
	printf("This was event number %i\n",i++);
	printf("==============================================================\n");
	if(maxEvents!=0 && i>maxEvents) return 0;
      }
  }
  
  return 0;
}
