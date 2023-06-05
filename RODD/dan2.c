#include "dan2.h"

void startup(void) {
    setupHbook(); 
    if (logfile != NULL) 
    logfl = fopen(logfile, "a");
  else
    logfl = fopen(default_log, "a");
  now = time(NULL); 
}

static void 
shutdown() {
  now = time(NULL);
  closeHbook();  
  fprintf(stdout,"dan2 processed %d events.\n\n",number_of_events);
}

void
Usage(FILE *outfile){
  fprintf(outfile,"\n");
  fprintf(outfile,"usage of dan2 is as follows\n");
  fprintf(outfile,"\n");
}

void main(int nargs, char *args[]) {
  char *argptr;
  int i;
  fprintf(stdout,"dan2 data analysis successfully launched.\n");
  for (i = 1; i < nargs; i++) {
    if (*args[i] == '-') {
      argptr = args[i] + 1;
      switch(*argptr) {
      case 'D':
	debug_mode_engaged=1;
	fprintf(stdout,"Debugging mode engaged.\n");
	break;
      case 'h':
	Usage(stdout);
	break;
      default:
	fprintf(stderr, "Unknown argument : %s\n", argptr);
	Usage(stderr);
      } /* end switch */
    }
    else{
      if(n_data_files>=MAX_NUMBER_OF_DATA_FILES){
	fprintf(stderr,"I'm sorry, but the maximum number of data files is %d.\n",MAX_NUMBER_OF_DATA_FILES);
	exit(1);
      }
      else{ /* will not exceed the max number of input files */
	if((data_files[n_data_files]=fopen(args[i],"r"))==NULL){
	  fprintf(stderr,"Trouble opening file >%s<.\n",args[i]);
	}
	else{
	  fprintf(stdout,"Opened >%s< for input.\n",args[i]);
	  
	}
	n_data_files++;
      }
    }
  }
  fprintf(stdout,"Command line parsing finished.\n");
  if(n_data_files==1)
    fprintf(stdout,"Analyzing %d file.\n",n_data_files);
  else
    fprintf(stdout,"Analyzing %d files.\n",n_data_files);
  startup();
  fprintf(stdout,"Finished setup.\n");
  fflush(stdout);
  /* the main driving loop */
  analyze_files();

  shutdown();
}

void analyze_files(void){
  int iter;
  for(iter=0;iter<n_data_files;iter++){
    analyze_file(data_files[iter]);
  }
}

void analyze_file(FILE *in_file){
  int return_value,current_run;
  itape_header_t *current_event=NULL;
  
  /*  lgdSetup(3,0.2,0.05,LGDCLUS_LOG2,0); */
  
  current_event=(itape_header_t *)malloc(BUFSIZE);
  return_value=data_read(fileno(in_file),current_event,BUFSIZE);
  while((max_number_of_events ? number_of_events<max_number_of_events : 1)&&((return_value=data_read(fileno(in_file),current_event,BUFSIZE))!=DATAIO_EOF)){
    switch(return_value){      
    case DATAIO_OK:
      /* Got an event */
      /*      printf("Event number %d.\n",current_event->eventNo);*/
      switch(current_event->type){
      case TYPE_TAPEHEADER:
	/* Get the tape number */
	break;
      case TYPE_ITAPE:
	if(current_event->runNo != current_run){
	  /* Load databases */
	  if(configureIt(current_event)){
	    fprintf(stderr,"Configuration failed!\n");
	    exit(5);
	  }
	  else{
	    printf("configureIt suceeded.\n");
	  }
	  current_run = current_event->runNo;
	}
	do_event(current_event);
	break;
      default:
	/* Oh no, no idea what this is.... */
	fprintf(stderr,"doAnEvent: got event of unknown type : %d\n",current_event->type);
	break;
      }
      break;
    case DATAIO_EOF:
      /* EOF marker on tape */
      fprintf(stderr,"doAnEvent: EOF on tape\n");
      printf("End of file marker.\n");
      return;
    case DATAIO_EOT:
      /* Real EOF tape or file */
      printf("End of Tape marker.\n");
      return;
    case DATAIO_ERROR:
      fprintf(stderr,"doAnEvent: data_read() returned DATAIO_ERROR, Skipping event\n");
    case DATAIO_BADCRC:
      fprintf(stderr,"doAnEvent: data_read() returned DATAIO_BADCRC, Skipping event\n");
      
    }
    consider_event_number(stdout,number_of_events,1000);
  }
}

void consider_event_number(FILE *print_file, int N, int divisor){
  if(!(N%divisor)){
    fprintf(print_file,"Event %d.\n",N);
  }
}


int configureIt(itape_header_t *event){
  int actualRun;
  int index;
    
  if(umap_loadDb(event->runNo))
    return(1);
  if(lgdGeomSetup(event->runNo))
    return(1);
  
  lgdGetGeomParam(&nChannels,&nRows,&nCols);
  
  /* Load pedestals */
  
  if(ped)
    free(ped);
  ped = malloc(nChannels*sizeof(float));
  if(threshold)
    free(threshold);
  threshold = malloc(nChannels*sizeof(float));
  
  if(pedestalLoad(PED_LGD,event->runNo,nChannels,ped,threshold,&actualRun))
    return(1);
  
  if(!cc){
    cc = malloc(nChannels*sizeof(float));
    
    if(ccLoad(CC_LGD,event->runNo,nChannels,cc,&actualRun)){
      if(debug_mode_engaged){
	fprintf(stdout,"Help!  There was a problem loading calibration constants!\n\n");
	exit(8);
      }
      /*    if(ccFile){
	    FILE *fp;
	    
	    fprintf(stderr,"Read cc's from %s\n",ccFile);
	    fp = fopen(ccFile,"r");
	    for(index=0;index<nChannels;index++){
	    fscanf(fp,"%f",&cc[index]);
	    }
	    fclose(fp);
	    
	    } replace with map-cc getting code 24SEP97 CPS*/
      fprintf(stderr,"Set cc's to %f\n",baseCC);
      for(index=0;index<nChannels;index++){
	cc[index] = baseCC;
      }
    }
  }
  /* Setup the clusterizer */
  
  /*  lgdSetup(3,0.2,0.05,LGDCLUS_LOG2,0); */
  lgdSetup(3,0.35,0.15,LGDCLUS_LOG2,0);
  printf("\nFinished lgdSetup.  Whoohoo!\n\n");  
  return(0);
}

int do_event(itape_header_t *c_event){
  int iter1;
  adc_values_t *rpd_adcs=NULL;
  if(debug_mode_engaged){
    printf("\n\ndo_event on event number %d.\n",c_event->eventNo); 
  }
  if( (c_event->eventType == EV_DATA) && (c_event->trigger == TRIG_DATA)){
    if(debug_mode_engaged){
      fprintf(stdout,"Event is data, unpacking.\n");
    }
    if(unpackEvent(c_event,BUFSIZE))
      return 0;
    number_of_events++;
    
    if(!(rpd_adcs=data_getGroup(c_event,GROUP_RPD_ADCS,NULL))){
      printf("Failing to read rpd_adcs group: aborting\n\n");
      exit(100);
    }
    for(iter1=0;iter1<rpd_adcs->nadc;iter1++){
      hf1(11000+rpd_adcs->adc[iter1].channel,
	  rpd_adcs->adc[iter1].value,
	  1.0);
    }
    
    
    
  } /* if( (c_event->eventType == EV_DATA) &&
       (c_event->trigger == TRIG_DATA)) */
  
  if(debug_mode_engaged){
    fprintf(stdout,"event finished.\n");
  }
}

void closeHbook(void){
  hrput(0,histo_file_name,"N");
}

void setupHbook(void){
  int iter;
  char name_string[80];
  hlimit(HBOOK); 
  hbook1(1,"trigger bit set, all events",64,-0.5,63.5,0.0);
  
  for(iter=0;iter<N_RPD_CHANNELS;iter++){
    sprintf(name_string,"ADC, RPD channel %d",iter);
    hbook1(11000+iter,name_string,450,0.0,4095.0,0.0);
  }
}
