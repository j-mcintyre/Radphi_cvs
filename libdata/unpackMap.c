/*
  $Log: unpackMap.c,v $
  Revision 1.10  1999/06/17 17:33:58  radphi
  Changes made by radphi@jlabs1
    BSD now has long and short-range TDC groups; need to allow
      umapGetGroup to know about moduel models
     D.S. Armstrong

  Revision 1.9  1999/06/15 20:56:39  radphi
     added BSD and BGV detectors to hardware map and data groups
     D.S. Armstrong

  Revision 1.8  1998/06/26 23:56:46  radphi
  Changes made by radphi@jlabs1
  added LRS 1875a decoding  D.S. Armstrong

  Revision 1.7  1998/04/01 17:37:26  radphi
  Changes made by jonesrt@jlabs3
  bringing the code up to date with the database - R.T.Jones

 * Revision 1.6  1998/03/31  02:35:52  radphi
 * Changes made by jonesrt@zeus
 * new types for new readout modules - R.T.Jones
 *
 * Revision 1.5  1997/05/20  09:05:49  radphi
 * Changes made by lfcrob@jlabs2
 * Put in changes to properly handle scalers
 *
 * Revision 1.4  1997/05/20  04:36:17  radphi
 * Changes made by radphi@jlabs2
 * Added CPV TDC's, ESUM ADC's
 *
 * Revision 1.3  1997/05/14  05:43:26  radphi
 * Changes made by lfcrob@jlabs2
 * Added routine to extract HW info based on
 * logical channel number
 *
 * Revision 1.2  1997/04/06  14:54:20  lfcrob
 * Fixed group for tagger TDC's
 *
 * Revision 1.1  1997/02/26  19:59:10  lfcrob
 * Initial revision
 *
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <disData.h>
#include <param.h>
#include <umap.h>


static char rcsid[] = "$Id: unpackMap.c,v 1.10 1999/06/17 17:33:58 radphi Exp $";

static void freeUmap(void);
static int umapMapModel(const char *model);
static int umapMapType(const char *type);
static int umapMapBus(const char *bus);
static int umapMapDetector(const char *detector);
static int umapGetGroup(int detector, int type, int model);

static umap_t *umap=NULL;

int umap_getHW(int detector, int module, int channel,int *slot,int *hwChannel)
{
  int rocNo;
  int slotNo;
  int channelNo;

  for(rocNo=0;rocNo<umap->nRocs;rocNo++)
    for(slotNo=0;slotNo<umap->roc[rocNo]->nModules;slotNo++)
      if(umap->roc[rocNo]->module[slotNo]->moduleType == module){
	for(channelNo=0;channelNo <umap->roc[rocNo]->module[slotNo]->nChannels;channelNo++){
	  if((umap->roc[rocNo]->module[slotNo]->channel[channelNo].detector == detector) &&
	     (umap->roc[rocNo]->module[slotNo]->channel[channelNo].swChannel == channel) ){
	    *slot = umap->roc[rocNo]->module[slotNo]->slot;
	    *hwChannel = channelNo;
	    return(0);
	  }	  
	}
      }
  return(1);
}

int umap_getInfo(int roc, int slot, int channel, umapData_t *data)
{
  int rocNo;
  int slotNo;
  int channelNo;

  for(rocNo=0;rocNo<umap->nRocs;rocNo++)
    if(umap->roc[rocNo]->rocNum == roc){
      for(slotNo=0;slotNo<umap->roc[rocNo]->nModules;slotNo++)
	if(umap->roc[rocNo]->module[slotNo]->slot == slot){
	  data->swChannel = umap->roc[rocNo]->module[slotNo]->channel[channel].swChannel;
	  data->group = umap->roc[rocNo]->module[slotNo]->channel[channel].group;
	  data->model = umap->roc[rocNo]->module[slotNo]->model;
	  data->type = umap->roc[rocNo]->module[slotNo]->moduleType;
	  data->detector = umap->roc[rocNo]->module[slotNo]->channel[channel].detector;
	  data->bus= umap->roc[rocNo]->module[slotNo]->bus;
	  return(0);
	}
    }
  return(1);
}

int umap_loadDb(int runNo)
{
  int nRocs;
  int roc;
  char dbFile[1024];
  char tmp[1024];

  if(param_loadDatabase("umap-master",runNo)){
    fprintf(stderr,"umap_loadDb : Cannot load master database!\n");
    return(1);
  }
  nRocs=param_getInt("nrocs",0);
  if(umap){
    freeUmap();
  }
  umap=malloc(sizeof(umap_t));
  umap->nRocs=nRocs;
  umap->roc = malloc(sizeof(umapRoc_t)*nRocs);
  for(roc=0;roc<nRocs;roc++){
    sprintf(dbFile,"umap-roc%d",roc);
    if(param_loadDatabase(dbFile,runNo)){
      fprintf(stderr,"umap_loadDb : Cannot load database %s!\n",dbFile);
      return(1);
    }
    /* We have the ROC db */
    {
      int nModules;
      int rocNo;
      const char *type;
      const char *model;
      int module;
      const char *bus;
      int nChannels;
      int channel;
      sprintf(tmp,"roc.%d.nmodules",roc);
      nModules=param_getInt(tmp,0);
      sprintf(tmp,"roc.%d.address",roc);
      rocNo=param_getInt(tmp,0);
      
      /* Now we can malloc the roc */

      umap->roc[roc] = malloc(sizeof(umapRoc_t));
      umap->roc[roc]->nModules = nModules;
      umap->roc[roc]->rocNum = rocNo;
      umap->roc[roc]->module=malloc(sizeof(umapModule_t *)*nModules);

      for(module=0;module<nModules;module++){
	int dChannel;

	sprintf(tmp,"roc.%d.module.%d.type",roc,module);
	model = param_getValue(tmp);
	sprintf(tmp,"%s.bus",model);
	bus = param_getValue(tmp);
	sprintf(tmp,"%s.type",model);
	type = param_getValue(tmp);
	sprintf(tmp,"%s.nchannels",model);
	nChannels=param_getInt(tmp,0);

	umap->roc[roc]->module[module] = malloc(sizeof(umapModule_t));
	umap->roc[roc]->module[module]->nChannels = nChannels;
	umap->roc[roc]->module[module]->channel=malloc(sizeof(umapChannel_t)*nChannels);
	umap->roc[roc]->module[module]->model=umapMapModel(model);
	umap->roc[roc]->module[module]->moduleType=umapMapType(type);
	umap->roc[roc]->module[module]->bus=umapMapBus(bus);
	sprintf(tmp,"roc.%d.module.%d.slot",roc,module);
	umap->roc[roc]->module[module]->slot = param_getInt(tmp,0);

	/*	printf("Module %d: %s -  %d Channel %s %s\n",module,model,nChannels,bus,type);*/
	for(channel=0;channel<nChannels;channel++){
	  sprintf(tmp,"roc.%d.module.%d.channel.%d.detector",roc,module,channel);
	  /*	  detector = param_getValue(tmp);*/
	  umap->roc[roc]->module[module]->channel[channel].detector = umapMapDetector(param_getValue(tmp));
	  sprintf(tmp,"roc.%d.module.%d.channel.%d.channel",roc,module,channel);
	  umap->roc[roc]->module[module]->channel[channel].swChannel = param_getInt(tmp,0);
	  umap->roc[roc]->module[module]->channel[channel].group=umapGetGroup(umap->roc[roc]->module[module]->channel[channel].detector,umap->roc[roc]->module[module]->moduleType,umap->roc[roc]->module[module]->model);
	  
	}
      }
    }
    
    
  }
  return(0);
}

static int umapMapModel(const char *model)
{
 if(!strcmp(model,"LRS_1885M"))
   return(UMAP_LRS_1885M);
 else if(!strcmp(model,"STK_10C6"))
   return(UMAP_STK_10C6);
 else if(!strcmp(model,"LRS_1151"))
   return(UMAP_LRS_1151);
 else if(!strcmp(model,"STK_7200"))
   return(UMAP_STK_7200);
 else if(!strcmp(model,"IU_01e0"))
   return(UMAP_IU_ADC);
 else if(!strcmp(model,"IU_ADC"))
   return(UMAP_IU_ADC);
 else if(!strcmp(model,"IU_MAM"))
   return(UMAP_IU_MAM);
 else if(!strcmp(model,"LRS_1877"))
   return(UMAP_LRS_1877);
 else if(!strcmp(model,"LRS_1875"))
   return(UMAP_LRS_1875);

 fprintf(stderr,"umapMapModel: Unknown module: %s\n",model);
 return(UMAP_UNKNOWN);
}

static int umapMapDetector(const char *detector)
{
 if(!strcmp(detector,"LGD"))
   return(UMAP_LGD);
 else if(!strcmp(detector,"ESUM"))
   return(UMAP_ESUM);
 else if(!strcmp(detector,"RPD"))
   return(UMAP_RPD);
 else if(!strcmp(detector,"CPV"))
   return(UMAP_CPV);
 else if(!strcmp(detector,"UPV"))
   return(UMAP_UPV);
 else if(!strcmp(detector,"SCALER"))
   return(UMAP_SCALER);
 else if(!strcmp(detector,"TAGGER"))
   return(UMAP_TAGGER);
 else if(!strcmp(detector,"BSD"))
   return(UMAP_BSD);
 else if(!strcmp(detector,"BGV"))
   return(UMAP_BGV);
 else if(!strcmp(detector,"UNUSED"))
   return(UMAP_UNUSED);
 fprintf(stderr,"umapMapDetector: Unknown detector: %s\n",detector);
 return(UMAP_UNKNOWN);
}

static int umapMapBus(const char *bus)
{
 if(!strcmp(bus,"CAMAC"))
   return(UMAP_CAMAC);
 else if(!strcmp(bus,"VME"))
   return(UMAP_VME);
 else if(!strcmp(bus,"FASTBUS"))
   return(UMAP_FASTBUS);
 fprintf(stderr,"umapMapBus: Unknown bus: %s\n",bus);
 return(UMAP_UNKNOWN);
}
static int umapMapType(const char *type)
{
 if(!strcmp(type,"ADC"))
   return(UMAP_ADC);
 else if(!strcmp(type,"TDC"))
   return(UMAP_TDC);
 else if(!strcmp(type,"SCALER"))
   return(UMAP_SCALER);
 else if(!strcmp(type,"REG"))
   return(UMAP_REG);
 fprintf(stderr,"umapMapType: Unknown module type: %s\n",type);
 return(UMAP_UNKNOWN);
}
static int umapGetGroup(int detector, int type, int model)
{
  switch(detector){
  case UMAP_LGD:
    switch(type){
    case UMAP_ADC:
      return(GROUP_LGD_ADCS);
    case UMAP_REG:
      return(GROUP_MAM_REGS);
    case UMAP_TDC:
      fprintf(stderr,"umapGetGroup : Don't know group for LGD TDCs!\n");
      return(GROUP_NOGROUP);
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for LGD scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  case UMAP_ESUM:
    switch(type){
    case UMAP_ADC:
      return(GROUP_ESUM_ADCS);
    case UMAP_TDC:
      fprintf(stderr,"umapGetGroup : Don't know group for ESUM TDCs!\n");
      return(GROUP_NOGROUP);
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for ESUM scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
    
  case UMAP_RPD:
    switch(type){
    case UMAP_ADC:
      return(GROUP_RPD_ADCS);
    case UMAP_TDC:
      return(GROUP_RPD_TDCS);
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for RPD scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  case UMAP_TAGGER:
    switch(type){      
    case UMAP_ADC:
      fprintf(stderr,"umapGetGroup : Don't know group for Tagger ADCs!\n");
      return(GROUP_NOGROUP);
    case UMAP_TDC:
      return(GROUP_TAGGER_TDCS);
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for Tagger scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  case UMAP_CPV:
    switch(type){
    case UMAP_ADC:
      return(GROUP_CPV_ADCS);
    case UMAP_TDC:
      return(GROUP_CPV_TDCS);
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for CPV scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  case UMAP_UPV:
    switch(type){
    case UMAP_ADC:
      return(GROUP_UPV_ADCS);
    case UMAP_TDC:
      return(GROUP_UPV_TDCS);
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for UPV scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  case UMAP_BSD:
    switch(type){
    case UMAP_ADC:
      return(GROUP_BSD_ADCS);
    case UMAP_TDC:
      switch(model){
      case UMAP_LRS_1877:
	return(GROUP_BSD_TDCS_LONG);
      case UMAP_LRS_1875:
	return(GROUP_BSD_TDCS_SHORT);
      default:
	fprintf(stderr,"umapGetGroup : wrong TDC module for BSD: %d\n",model);
      }
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for BSD scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  case UMAP_BGV:
    switch(type){
    case UMAP_ADC:
      return(GROUP_BGV_ADCS);
    case UMAP_TDC:
      return(GROUP_BGV_TDCS);
    case UMAP_SCALER:
      fprintf(stderr,"umapGetGroup : Don't know group for BGV scalers!\n");
      return(GROUP_NOGROUP);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  case UMAP_SCALER:
    switch(type){
    case UMAP_ADC:
      fprintf(stderr,"umapGetGroup : Don't know group for scaler ADCs!\n");
      return(GROUP_NOGROUP);
    case UMAP_TDC:
      fprintf(stderr,"umapGetGroup : Don't know group for scaler TDCs!\n");
      return(GROUP_NOGROUP);
    case UMAP_SCALER:
      return(GROUP_SCALERS);
    default:
      fprintf(stderr,"umapGetGroup : Unknown module type: %d\n",type);
      return(GROUP_NOGROUP);
    }
  default:
    break;
  }
  return(GROUP_NOGROUP);
}
static void freeUmap(void)
{
  int roc;
  int module;
  int channel;

  /* Gotta go backwards....<sigh> */

  for(roc=0;roc<umap->nRocs;roc++){
    for(module=0;module<umap->roc[roc]->nModules;module++){
      free(umap->roc[roc]->module[module]->channel);
      free(umap->roc[roc]->module[module]);
    }
    free(umap->roc[roc]->module);
    free(umap->roc[roc]);
  }
  free(umap->roc);
  free(umap);
  umap=NULL;
  return;
}

void umapDataPrint(FILE *fp,umapData_t *data)
{
  char *model=NULL;
  char *type=NULL;
  char *detector=NULL;
  char *group=NULL;
  char *bus;

  switch(data->group){
  case GROUP_LGD_ADCS:
    group = "GROUP_LGD_ADCS";
    break;
  case GROUP_RPD_ADCS:
    group = "GROUP_RPD_ADCS";
    break;
  case GROUP_RPD_TDCS:
    group = "GROUP_RPD_TDCS";
    break;
  case GROUP_CPV_ADCS:
    group = "GROUP_CPV_ADCS";
    break;
  case GROUP_CPV_TDCS:
    group = "GROUP_CPV_TDCS";
    break;
  case GROUP_UPV_ADCS:
    group = "GROUP_UPV_ADCS";
    break;
  case GROUP_UPV_TDCS:
    group = "GROUP_UPV_TDCS";
    break;
  case GROUP_BSD_ADCS:
    group = "GROUP_BSD_ADCS";
    break;
  case GROUP_BSD_TDCS_LONG:
    group = "GROUP_BSD_TDCS_LONG";
    break;
  case GROUP_BSD_TDCS_SHORT:
    group = "GROUP_BSD_TDCS_SHORT";
    break;
  case GROUP_BGV_ADCS:
    group = "GROUP_BGV_ADCS";
    break;
  case GROUP_BGV_TDCS:
    group = "GROUP_BGV_TDCS";
    break;
  case GROUP_TAGGER_TDCS:
    group = "GROUP_TAGGER_TDCS";
    break;
  case GROUP_SCALERS:
    group = "GROUP_SCALERS";
    break;
  default:
    group = "UNKNOWN";
    break;
  }
  switch(data->model){
  case UMAP_LRS_1885M:
    model = "LRS_1885M";
    break;
  case UMAP_STK_10C6:
    model = "STK_10C6";
    break;
  case UMAP_LRS_1151:
    model = "LRS_1151";
    break;
  case UMAP_STK_7200:
    model = "STK_7200";
    break;
  }
  switch(data->type){
  case UMAP_ADC:
    type = "ADC";
    break;
  case UMAP_TDC:
    type = "TDC";
    break;
  case UMAP_SCALER:
    type = "SCALER";
    break;
  }
  switch(data->bus){
  case UMAP_VME:
    bus = "VME";
    break;
  case UMAP_FASTBUS:
    bus = "FASTBUS";
    break;
  case UMAP_CAMAC:
    bus = "CAMAC";
    break;
  }
  switch(data->detector){
  case UMAP_LGD:
    detector = "LGD";
    break;
  case UMAP_RPD:
    detector = "RPD";
    break;
  case UMAP_SCALER:
    detector = "SCALER";
    break;
  case UMAP_CPV:
    detector = "CPV";
    break;
  case UMAP_UPV:
    detector = "UPV";
    break;
  case UMAP_BGV:
    detector = "BGV";
    break;
  case UMAP_BSD:
    detector = "BSD";
    break;
  case UMAP_TAGGER:
    detector = "TAGGER";
    break;
  case UMAP_UNUSED:
    detector = "N/A";
    break;
  }
  fprintf(fp,"%s (%s %s) for %s channel %d, stored in %s\n",model,bus,type,detector,data->swChannel,group);
}
