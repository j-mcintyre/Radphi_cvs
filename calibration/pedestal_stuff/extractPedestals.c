/*
  $Log: extractPedestals.c,v $
  Revision 1.7  2000/06/09 01:04:30  radphi
  Changes made by radphi@urs3
  added some debugging features

  Revision 1.6  2000/05/16 16:25:36  radphi
  Changes made by radphi@urs1
  now takes care of correct number of ADCs in the upper crate

  Revision 1.5  1999/07/29 10:40:53  radphi
  Changes made by gevmage@jlabs1

  Revision 1.4  1999/07/04 17:33:50  radphi
  Changes made by gevmage@jlabs1
  clarifying file name conventions

  Revision 1.3  1999/07/04 14:43:21  radphi
  Changes made by gevmage@jlabs1
  all 1999 detectors now have pedestals set up.  Also, the
  zero supression thresholds for the DAQ are being set.

  Revision 1.2  1999/07/03 16:49:26  radphi
  Changes made by gevmage@jlabs1
  updating, about to add 1999 detectors

  Revision 1.1.1.1  1998/06/09 01:02:20  radphi
  Initial revision by radphi@urs2

 * Revision 1.2  1997/07/28  11:52:26  radphi
 * Changes made by lfcrob@dustbunny
 * Added other remaining detectors
 *
 * Revision 1.1.1.1  1997/05/24  00:46:54  radphi
 * Initial revision by lfcrob@jlabs2
 *
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <disData.h>
#include <dataIO.h>
#include <unpackData.h>
#include <triggerType.h>
#include <eventType.h>
#include <itypes.h>
#include <lgdGeom.h>
#include <rstats.h>
#include <umap.h>

#define BUFSIZE 100000

#define TEMP_NUMBER_OF_LGD_CHANNELS 784
#define TEMP_NUMBER_OF_CPV_CHANNELS 30
#define TEMP_NUMBER_OF_BSD_CHANNELS 48
#define TEMP_NUMBER_OF_BGV_CHANNELS 48
#define TEMP_NUMBER_OF_UPV_CHANNELS 8
#define TEMP_NUMBER_OF_AUX_CHANNELS 10
#define ADCS_IN_TOP_CRATE 15
#define ADCS_IN_BOTTOM_CRATE 13
#define CHANNELS_PER_ADC 32
#define ADC_BLOCK_HEADER_WORD 0xda000adc
#define ADC_DUMMY_WORD 0xdaffffff
#define ADC_BLOCK_END_WORD 0xda000cda
#define ADC_DATA_MASK 0x00000fff

#define OS_LGD 0
#define OS_UPV 1000
#define OS_BSD 2000
#define OS_BGV 2200
#define OS_CPV 2400
#define OS_AUX 2600
#define OS_TOP 4000
#define OS_BOTTOM 5000

static const char rcsid[]="$Id: extractPedestals.c,v 1.7 2000/06/09 01:04:30 radphi Exp $";

void main(int argc, char *argv[]);
int processEvent(itape_header_t *event);


static int nEvents=0;
static int nPedestals=0;
static int *hasData=NULL;
static int nChannels,nRows,nCols;
int debug_mode_engaged=0;


void main(int argc, char *argv[])
{
  int iarg;
  char *argptr;
  char *inputFile=NULL;
  FILE *fp;
  int fn;
  int ret;
  itape_header_t *event=NULL;
  int running=1;
  int currentRun=0;
  int index;
  average_t av;
  int nHaveData=0;
  int threshold;
  int N_top,N_bottom;
  char fname[1024];
  FILE *fpOut;
  
  N_top=ADCS_IN_TOP_CRATE*CHANNELS_PER_ADC;
  N_bottom=ADCS_IN_BOTTOM_CRATE*CHANNELS_PER_ADC;
  
  for(iarg=1;iarg<argc;iarg++)
    if(*argv[iarg] == '-'){
      argptr = argv[iarg]+1;
      switch(*argptr){
      default:
	;
      }
    }
    else{
      inputFile = argv[iarg];
    }
  
  /* Open the input file */

  if(inputFile){
    fn=0;
    fn = evt_open(inputFile,"r");
    if(!fn){
      fprintf(stderr,"Failed to open %s for input\n",inputFile);
      exit(1);
    }
    fprintf(stderr,"Reading data from %s\n",inputFile);
  }
  event = malloc(BUFSIZE);

  while(running){
    ret=evt_data_read(fn,event,BUFSIZE);
    if(debug_mode_engaged)
      fprintf(stdout,"return value from evt_data_read is %d.\n",ret);
    fflush(stdout);
    switch(ret){
    case DATAIO_OK:
      if(event->runNo != currentRun){
	if(currentRun != 0){
	  fprintf(stderr,"Can only process one run at a time!\n");
	  break;
	}
	if(umap_loadDb(event->runNo)){	  
	  fprintf(stderr,"Failed to load umap database\n");
	  running=0;
	  break;
	}
	lgdGeomSetup(event->runNo);
	lgdGetGeomParam(&nChannels,&nRows,&nCols);
	if(hasData)
	  free(hasData);
	hasData = malloc(10000);
	memset(hasData,0,10000);
	for(index=0;index<nChannels;index++)
	  zero_av(index);
	currentRun=event->runNo;
      }
      if(debug_mode_engaged)
	fprintf(stdout,"Beginning processEvent.\n");
      processEvent(event);
      if( (nEvents%1000) == 0){
	printf("%d %d\r",nEvents,nPedestals);
	fflush(stdout);
      }
      break;
    case DATAIO_EOT:
    case DATAIO_EOF:
      fprintf(stdout,"At end of file.\n");
      running=0;
      break;
    default:
      fprintf(stderr,"Unkown return code from data_read: %d\n",ret);
      running=0;
      break;
    }
  }  

  /* ADCs top crate */

  sprintf(fname,"ped_adcs_top_%d_%d.dat",currentRun,N_top);
  if((fpOut = fopen(fname,"w")) == NULL){
    fprintf(stderr,"Failed to open %s",fname);
    exit(1);
  }
  fprintf(fpOut,"PED_ADCS_TOP %d %d\n",currentRun,N_top);
  for(index=0;index<N_top;index++){
    if(hasData[index+OS_TOP]){
      nHaveData++;
      get_av(index+OS_TOP,&av);
      if(av.sigma>5){
	fprintf(stdout,"top crate, slot %d, channel %d has a pedestal sigma of %5.2f!\n",
		(24-(index/32)),index%32,av.sigma);
      }
      /*      threshold = (int) (av.average+3*av.sigma); */
      threshold = (int) (av.average+5.0);
      fprintf(fpOut,"%d %f %d\n",index,av.average,threshold);
      if(av.average>35.0){
	fprintf(stdout,"top crate channel %d has pedestal %.2f\n",
		index,av.average);
	fprintf(stdout,"that is slot %d channel %d\n",
		(24-(index/32)),index%32);
      }
    }
  }
  fclose(fpOut);

  /* ADCs bottom crate */

  sprintf(fname,"ped_adcs_bottom_%d_%d.dat",currentRun,N_bottom);
  if((fpOut = fopen(fname,"w")) == NULL){
    fprintf(stderr,"Failed to open %s",fname);
    exit(1);
  }
  fprintf(fpOut,"PED_ADCS_BOTTOM %d %d\n",currentRun,N_bottom);
  for(index=0;index<N_bottom;index++){
    if(hasData[index+OS_BOTTOM]){
      nHaveData++;
      get_av(index+OS_BOTTOM,&av);
      if(av.sigma>5){
	fprintf(stdout,"bottom crate, slot %d, channel %d has a pedestal sigma of %5.2f !\n",
		(24-(index/32)),index%32,av.sigma);
      }
      /*      threshold = (int) (av.average+3*av.sigma);*/
      threshold = (int) (av.average+5.0);
      fprintf(fpOut,"%d %f %d\n",index,av.average,threshold);
      if(av.average>35.0){
	fprintf(stdout,"bottom crate channel %d has pedestal %.2f\n",
		index,av.average);
	fprintf(stdout,"that is slot %d channel %d\n",
		(20-(index/32)),index%32);
      }
    }
  }
  fclose(fpOut);

  /* LGD */

  sprintf(fname,"ped_lgd_%d_%d.dat",currentRun,nChannels);
  if((fpOut = fopen(fname,"w")) == NULL){
    fprintf(stderr,"Failed to open %s",fname);
    exit(1);
  }
  fprintf(fpOut,"PED_LGD %d %d\n",currentRun,nChannels);
  for(index=0;index<nChannels;index++){
    if(hasData[index+OS_LGD]){
      nHaveData++;
      get_av(index+OS_LGD,&av);
      /*      threshold = (int) (av.average+3*av.sigma); */
      threshold = (int) (av.average+5.0);
      fprintf(fpOut,"%d %f %d\n",index,av.average,threshold);
      if(av.average>35.0){
	fprintf(stdout,"LGD channel %d has pedestal %.2f\n",
		index,av.average);
	fprintf(stdout,"That is row %d column %d\n",
		index/28,index%28);
	fprintf(stdout,"ADC channel ID is %08x\n",
		((index/28)<<24)|((index%28)<<16));
      }
    }
  }
  fclose(fpOut);

  /* UPV */

  sprintf(fname,"ped_upv_%d_%d.dat",currentRun,TEMP_NUMBER_OF_UPV_CHANNELS);
  if((fpOut = fopen(fname,"w")) == NULL){
    fprintf(stderr,"Failed to open %s",fname);
    exit(1);
  }
  fprintf(fpOut,"PED_UPV %d %d\n",currentRun,TEMP_NUMBER_OF_UPV_CHANNELS);
  for(index=0;index<TEMP_NUMBER_OF_UPV_CHANNELS;index++){
    if(hasData[index+OS_UPV]){
      nHaveData++;
      get_av(index+OS_UPV,&av);
      /*      threshold = (int) (av.average+3*av.sigma);*/
      threshold = (int) (av.average+5.0); 
      fprintf(fpOut,"%d %f %d\n",index,av.average,threshold);
    }
  }
  fclose(fpOut);

  /* BSD */

  sprintf(fname,"ped_bsd_%d_%d.dat",currentRun,TEMP_NUMBER_OF_BSD_CHANNELS);
  if((fpOut = fopen(fname,"w")) == NULL){
    fprintf(stderr,"Failed to open %s",fname);
    exit(1);
  }
  fprintf(fpOut,"PED_BSD %d %d\n",currentRun,TEMP_NUMBER_OF_BSD_CHANNELS);
  for(index=0;index<TEMP_NUMBER_OF_BSD_CHANNELS;index++){
    if(hasData[index+OS_BSD]){
      nHaveData++;
      get_av(index+OS_BSD,&av);
      /*      threshold = (int) (av.average+3*av.sigma); */
      threshold = (int) (av.average+5.0);
      fprintf(fpOut,"%d %f %d\n",index,av.average,threshold);
    }
  }
  fclose(fpOut);

  /* BGV */

  sprintf(fname,"ped_bgv_%d_%d.dat",currentRun,TEMP_NUMBER_OF_BGV_CHANNELS);
  if((fpOut = fopen(fname,"w")) == NULL){
    fprintf(stderr,"Failed to open %s",fname);
    exit(1);
  }
  fprintf(fpOut,"PED_BGV %d %d\n",currentRun,TEMP_NUMBER_OF_BGV_CHANNELS);
  for(index=0;index<TEMP_NUMBER_OF_BGV_CHANNELS;index++){
    if(hasData[index+OS_BGV]){
      nHaveData++;
      get_av(index+OS_BGV,&av);
      /*      threshold = (int) (av.average+3*av.sigma); */
      threshold = (int) (av.average+5.0);
      fprintf(fpOut,"%d %f %d\n",index,av.average,threshold);
    }
  }
  fclose(fpOut);

  /* CPV */

  sprintf(fname,"ped_cpv_%d_%d.dat",currentRun,TEMP_NUMBER_OF_CPV_CHANNELS);
  if((fpOut = fopen(fname,"w")) == NULL){
    fprintf(stderr,"Failed to open %s",fname);
    exit(1);
  }
  fprintf(fpOut,"PED_CPV %d %d\n",currentRun,TEMP_NUMBER_OF_CPV_CHANNELS);
  for(index=0;index<TEMP_NUMBER_OF_CPV_CHANNELS;index++){
    if(hasData[index+OS_CPV]){
      nHaveData++;
      get_av(index+OS_CPV,&av);
      /*      threshold = (int) (av.average+3*av.sigma); */
      threshold = (int) (av.average+5.0);
      fprintf(fpOut,"%d %f %d\n",index,av.average,threshold);
    }
  }
  fclose(fpOut);
  printf("Processed %d events, %d pedestal events\n",nEvents,nPedestals);
  printf("Calculated pedestals for %d channels\n",nHaveData);
}



int processEvent(itape_header_t *event)
{
  adc_values_t *lgd=NULL;
  adc_values_t *upv=NULL;
  adc_values_t *bsd=NULL;
  adc_values_t *bgv=NULL;
  adc_values_t *cpv=NULL;
  uint32 *raw=NULL;
  int raw_group_size;
  int current_raw_position=0;
  int n_top_crate,n_bottom_crate;
  int index;
  int adc_val;

  n_top_crate=ADCS_IN_TOP_CRATE*CHANNELS_PER_ADC;
  n_bottom_crate=ADCS_IN_BOTTOM_CRATE*CHANNELS_PER_ADC;  
  if(event->eventType == EV_DATA){
    nEvents++;
    if(debug_mode_engaged)
      fprintf(stdout,"trigger=%d.\n",event->trigger);
    
    if(event->trigger == TRIG_PED){
      if(unpackEvent(event,BUFSIZE)){
	fprintf(stderr,"Failed to unpack event\n");
	return(1);
      }
      if(debug_mode_engaged){
	fprintf(stderr,"unpacking finished.\n");
      }
      if(debug_mode_engaged){
	fprintf(stdout,"Finished unpacking event.\n");
      }
      if((raw = data_getGroup(event,GROUP_RAW,NULL)) == NULL){
	fprintf(stderr,"No RAW data group!\n");
	return(1);
      }
      raw_group_size=data_getGroupSize(event,GROUP_RAW,NULL);
      if((lgd = data_getGroup(event,GROUP_LGD_ADCS,NULL)) == NULL){
	fprintf(stderr,"No LGD ADC's!\n");
	return(1);
      }
      if(debug_mode_engaged)
	fprintf(stdout,"Got LGD info.\n");
      if((upv = data_getGroup(event,GROUP_UPV_ADCS,NULL)) == NULL){
	fprintf(stderr,"No UPV ADC's!\n");
	return(1);
      }
      if((bsd = data_getGroup(event,GROUP_BSD_ADCS,NULL)) == NULL){
	fprintf(stderr,"No BSD ADC's!\n");
	return(1);
      }
      if((bgv = data_getGroup(event,GROUP_BGV_ADCS,NULL)) == NULL){
	fprintf(stderr,"No BGV ADC's!\n");
	return(1);
      }
      if((cpv = data_getGroup(event,GROUP_CPV_ADCS,NULL)) == NULL){
	fprintf(stderr,"No CPV ADC's!\n");
	return(1);
      }
      nPedestals++;
      /* Now do something with it! */
      
      if(debug_mode_engaged)
	fprintf(stdout,"Beginning this event's accumulation.\n");
      if(find_next_adc_block(raw,&current_raw_position,raw_group_size)){
	fprintf(stderr,"There has been a problem finding first ADC block.\n");
	return(1);
      }
      for(index=0;index<n_top_crate;index++){
	hasData[index+OS_TOP]++;
	if(get_next_adc_value(raw,&current_raw_position,
			      raw_group_size,&adc_val)){
	  fprintf(stderr,"problem getting the next ADC value; top crate.\n");
	  return(1);
	}
	accumulate_av(index+OS_TOP,adc_val);
      }
      if(find_next_adc_block(raw,&current_raw_position,raw_group_size)){
	fprintf(stderr,"There has been a problem finding second adc block.\n");
	return(1);
      }
      for(index=0;index<n_bottom_crate;index++){
	hasData[index+OS_BOTTOM]++;
	if(get_next_adc_value(raw,&current_raw_position,
			      raw_group_size,&adc_val)){
	  fprintf(stderr,"problem getting the next ADC value; bottom crate.\n");
	  return(1);
	}
	accumulate_av(index+OS_BOTTOM,adc_val);
      }
      for(index=0;index<lgd->nadc;index++){
	hasData[lgd->adc[index].channel+OS_LGD]++;
	accumulate_av(lgd->adc[index].channel+OS_LGD,lgd->adc[index].value);
      }
      for(index=0;index<upv->nadc;index++){
	hasData[upv->adc[index].channel+OS_UPV]++;
	accumulate_av(upv->adc[index].channel+OS_UPV,upv->adc[index].value);
      }
      for(index=0;index<bsd->nadc;index++){
	hasData[bsd->adc[index].channel+OS_BSD]++;
	accumulate_av(bsd->adc[index].channel+OS_BSD,bsd->adc[index].value);
      }
      for(index=0;index<bgv->nadc;index++){
	hasData[bgv->adc[index].channel+OS_BGV]++;
	accumulate_av(bgv->adc[index].channel+OS_BGV,bgv->adc[index].value);
      }
      for(index=0;index<cpv->nadc;index++){
	hasData[cpv->adc[index].channel+OS_CPV]++;
	accumulate_av(cpv->adc[index].channel+OS_CPV,cpv->adc[index].value);
      }
      return(0);
    }
    if(debug_mode_engaged)
      fprintf(stdout,"Wasn't a pedstal event.\n");
    
  }
  return(0);
}

int get_next_adc_value(long unsigned *raw_data_buffer,
		       int *current_position, int data_length,
		       int *return_data){
  (*current_position)++;
  if((*current_position)>=data_length){
    fprintf(stderr,"get_next_adc_value ran off the end of the raw data.\n");
    return(1);
  }
  *return_data=( (*(raw_data_buffer+(*current_position)))&ADC_DATA_MASK);
  return(0);
}
  
int find_next_adc_block(long unsigned *raw_data_buffer,
			int *current_position, int data_length){
  int found_header=0;
  int i;
  if((*(raw_data_buffer+(*current_position)))==ADC_BLOCK_HEADER_WORD)
    found_header=1;
  while(!found_header){
    (*current_position)++;
    if((*current_position)>=data_length){
      fprintf(stderr,"find_next_adc_block ran off the end of the event.\n");
      return(1);
    }
    /*    fprintf(stdout,"checking word: %#010.8x\n",(*(raw_data_buffer+(*current_position))) ); */
    if((*(raw_data_buffer+(*current_position)))==ADC_BLOCK_HEADER_WORD)
      found_header=1;
  }
  (*current_position)++;
  if((*current_position)>=data_length){
    fprintf(stderr,"find_next_adc_block ran off the end of the event.\n");
    return(1);
  }
  if((*(raw_data_buffer+(*current_position)))==ADC_DUMMY_WORD){
    return(0);
  }
  else{
    fprintf(stderr,"%#010.8x is not the dummy word!\n",
	    ( *(raw_data_buffer+(*current_position)) ));
    fprintf(stderr,"Buffer miniblock:\n");
    for(i=-3;i<1;i++){
      fprintf(stderr,"\t%#010.8x\n",
	      (*(raw_data_buffer+(*current_position)+i)));
    }
    return(1);
  }
}
