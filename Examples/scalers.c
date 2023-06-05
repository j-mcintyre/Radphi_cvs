#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <disData.h>
#include <dataIO.h>
#include <itypes.h>
#include <param.h>
#include <eventType.h>
#include <triggerType.h>

#include "scalers.h"

static scatable *names;

scatable *get_scaler_names (int runNo)
{
  int i;
  int nRocs;
  int roc;
  char dbFile[1024];
  char recID[1024];
  static scatable *list=NULL;

  if(list!=NULL){
    for(i=0;i<MAX_SCALERS;i++)
      if(list->scaler[i]!=NULL)
        free(list->scaler[i]);
    free(list);
  }
  list = malloc(sizeof(scatable));
  for(i=0;i<MAX_SCALERS;i++)
    list->scaler[i] = NULL;
  if(param_loadDatabase("umap-master",runNo)) {
    fprintf(stderr,"get_scaler_names : Cannot load master database!\n");
    return NULL;
  }
  nRocs = param_getInt("nrocs",0);
  for(roc=0;roc<nRocs;roc++) {
    sprintf(dbFile,"umap-roc%d",roc);
    if(param_loadDatabase(dbFile,runNo)){
      fprintf(stderr,"get_scaler_names : Cannot load database %s!\n",dbFile);
      return NULL;
    }
    /* We have the ROC db */
    {
      int nModules;
      const char *type;
      const char *model;
      const char *bus;
      int module;
      int nChannels;
      int channel;
      sprintf(recID,"roc.%d.nmodules",roc);
      nModules=param_getInt(recID,0);

      for(module=0;module<nModules;module++){
	int dChannel;

	sprintf(recID,"roc.%d.module.%d.type",roc,module);
	model = param_getValue(recID);
	sprintf(recID,"%s.bus",model);
	bus = param_getValue(recID);
	sprintf(recID,"%s.type",model);
	type = param_getValue(recID);
	sprintf(recID,"%s.nchannels",model);
	nChannels=param_getInt(recID,0);

        if(!strcmp(type,"SCALER")) {
	  for(channel=0;channel<nChannels;channel++){
	    sprintf(recID,"roc.%d.module.%d.channel.%d.channel",
                           roc,module,channel);
	    if((i=param_getInt(recID,0)) >= MAX_SCALERS) {
              fprintf(stderr,"get_scaler_names : scaler channel %d",i);
              fprintf(stderr," in database greater than MAX_SCALERS!\n");
              return NULL;
            }
            if(list->scaler[i]!=NULL) {
              fprintf(stderr,"get_scaler_names : ");
              fprintf(stderr,"channel %d is duplicated in database!\n",i);
              return NULL;
            }
            list->scaler[i] = malloc(sizeof(scalabel));
	    sprintf(recID,"roc.%d.module.%d.channel.%d.name",
                           roc,module,channel);
            list->scaler[i]->tag = param_getValue(recID);
	    sprintf(recID,"roc.%d.module.%d.channel.%d.description",
                           roc,module,channel);
            list->scaler[i]->string = param_getValue(recID);
          }
	}
      }
    }
  }
  list->nscalers = 0;
  for(i=0;i<MAX_SCALERS;i++) {
    if(list->scaler[i]==NULL)
      list->nscalers = (list->nscalers==0) ? i : list->nscalers;
    else
      list->nscalers = (list->nscalers>0) ? -1 : list->nscalers;
  }
  if(list->nscalers < 0){
    fprintf(stderr,"get_scaler_names : ");
    fprintf(stderr,"the database contains gaps in the scaler channel list!\n");
    return NULL;
  }
  return list;
}

void make_scaler_panel (int ntid)
{
  FILE *fp;
  int i;

  fp = fopen("spanl.kumac","w");
  if (fp==NULL){
    perror("Could not open scratch file spanl.kumac");
    fprintf(stderr,"Please check the file system\n");
    return;
  }
  fprintf(fp,"*\n* Paw++ panel for looking at scalers\n*\n");
  fprintf(fp,"panel 0\n");

  /* Struck NIM */

  for (i=0;i<names->nscalers;i++)
    fprintf(fp,"panel %d.%2.2d 'exec scalers#lookat var=%s title=''%s''' '%s'\n",
            (i % 32)+1, (i / 32)+1,
            names->scaler[i]->tag,
            names->scaler[i]->string,
            names->scaler[i]->string);
  fprintf(fp,"panel 0 d 'Scalers' 341x339+869+585\n");
  fclose(fp);
}

int scalers_ (char *inspec, int *ntid, int *count, int *start, int len)
{
  int fd;
  int ret;
  int *raw;
  static int srun[3];
  static int *header=&srun[0];
  static int *scalArray=NULL;
  static itape_header_t *event=NULL;
  int records=0, srecords=0;
  scaler_values_t *scalers;
  int i;
  char infile[40];
  char *fin;
  
  if (!event) {
    event = malloc(BUFSIZE);
    if (scalArray == NULL)
      scalArray = malloc(sizeof(int)*MAX_SCALERS);
  }

  if (strncmp(inspec,"init",4) == 0) {
    struct rlimit rlim;
    getrlimit(RLIMIT_CORE,&rlim);
    rlim.rlim_cur = 1000;
    setrlimit(RLIMIT_CORE,&rlim);
    return -1;
  }

  strncpy(infile,inspec,len);
  if ((fin=(char *)index(infile,' ')) != NULL)
    fin[0]=0;
  fd = evt_open(infile,"r");
  if (fd==0)
    return -1;

  if (*start > 0) {
     printf("Skipping %d events...\n",*start);
     while  ( ((ret=evt_data_read(fd,event,BUFSIZE)) == DATAIO_OK)
            && (records++ < *start) ) {}
  }
  srun[0]=srun[1]=srun[2]=0;
  while  ( ((ret=evt_data_read(fd,event,BUFSIZE)) == DATAIO_OK)
         && (srecords < *count) ) {
    if ((event->eventType == EV_DATA) && (event->trigger != TRIG_DATA)) {
      if (event->runNo != srun[0]) {
        if (umap_loadDb(event->runNo)) {
          fprintf(stderr,"Warning: database access failure\n");
          continue;
        }
        if ((names=get_scaler_names(event->runNo)) == NULL) {
          fprintf(stderr,"scalers error: no scaler names found in database for run %d\n",
            event->runNo);
          return -1;
        }
        make_scaler_panel(*ntid);
        declare_scaler_ntuple(*ntid,names,header,scalArray);
        for (i=0;i<names->nscalers;i++)
          scalArray[i] = 0;
        srun[0] = event->runNo;
        srun[1] = event->eventNo;
        srun[2] = 0;
      }
      unpackEvent(event,BUFSIZE);
      if ((scalers=data_getGroup(event,GROUP_SCALERS,0)) != NULL) {
        if (scalers->nscalers != names->nscalers)
          fprintf(stderr,
                 "scalers warning: found %d scalers where %d expected!\n",
                  scalers->nscalers, names->nscalers);
        else
          for (i=0;i<scalers->nscalers;i++) {
            scalArray[scalers->scaler[i].channel] = 
                                         scalers->scaler[i].value;
	}
        srun[0] = event->runNo;
        srun[1] = event->eventNo;
        hfnt_(ntid);
        ++srun[2];
        ++srecords;
        {
          int c=srecords;
          float cf=srecords;
          while (c > 10) {
            cf /= 10;
            c /= 10;
          }
          if (cf == c) {
            printf("found scaler record %d of run %d after record %d\n",
                    srun[2],srun[0],records);
          }
        }
      }
    }
    ++records;
  }
  printf("\n%d data records processed, %d scaler records seen\n",
          records,srecords);
  evt_close(fd);
  return *count - srecords;
}
