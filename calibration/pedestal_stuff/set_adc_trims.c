#include <stdio.h>
#include <string.h>

#define ADC_PEDESTAL_SLOPE 0.891
#define TRIM_CHANGE_THROTTLE (1.0)
/* this is an attempt to throttle down the pedestal trim changes to prevent trims
   that won't settle down.  Change made June 6 2000 */

void usage(void);
int generate_new_trim(int target_value,
		      float current_pedestal,
		      int current_trim);
void usage(void){
  fprintf(stdout,"usage of set_adc_trims is as follows:\n");
  fprintf(stdout,"set_adc_trims <trim file> <ped file> <target ped>\n");
  fprintf(stdout,"<config file> is the pedestal trim config file name containing the \n");
  fprintf(stdout,"  old trims.  Typically, this is adc_trim_topcrate.config\n");
  fprintf(stdout,"  or some such.\n");
  fprintf(stdout,"<ped file> is the name of a file of produced by an extractPedestals call\n");
  fprintf(stdout,"  on data produced with <config file> in place.\n");
  fprintf(stdout,"  For example, ped_adcs_top_6022_480.dat.  It is important\n");
  fprintf(stdout,"  they are for the same crate (top or bottom) and the run\n");
  fprintf(stdout,"  that produced the pedestals used _that_ trim file.\n");
  fprintf(stdout,"  Otherwise, results are unpredictable.  \n");
  fprintf(stdout,"<target_ped> is an integer where you would like the new\n");
  fprintf(stdout,"  pedestals to be.\n");
}

void main(int nargs, char *args[]){
  int set_length,running_length=0;
  int i;
  char ID_string[50],run_string[6];
  FILE *trim_file=NULL;
  FILE *ped_file=NULL;
  FILE *output_file=NULL;
  char output_file_name[128];

  float target_value;
  int channel_counter=-1;
  int in_ped_num,in_thresh,in_slot,in_chan,in_trim;
  float in_ped;
  int out_trim;

  if(nargs!=4){
    usage();
    exit(1);
  }
  trim_file=fopen(args[1],"r");
  if(trim_file==NULL){
    fprintf(stderr,"There was a problem opening %s as trim file, exiting.\n",
	    args[1]);
    exit(1);
  } 
  ped_file=fopen(args[2],"r");
  if(ped_file==NULL){
    fprintf(stderr,"There was a problem opening %s as pedestal file, exiting.\n",
	    args[2]);
    exit(1);
  }
  fprintf(stdout,"Converting %s \n",args[3]);
  target_value=(float)atoi(args[3]);
  fprintf(stdout,"Files opened successfully.\n");
  fscanf(ped_file,"%s %s %d",ID_string,run_string,&set_length);
  if(!strcmp("PED_ADCS_TOP",ID_string))
    sprintf(output_file_name,"adc_trim_topcrate.config.new");
  else if(!strcmp("PED_ADCS_BOTTOM",ID_string))
    sprintf(output_file_name,"adc_trim_bottomcrate.config.new");
  else{
    fprintf(stdout,"I'm sorry, but your pedestal file does not have the\n");
    fprintf(stdout,"proper preamble.  The preamble I find is:\n");
    fprintf(stdout,"%s %s %d\n",ID_string,run_string,set_length);
    exit(1);
  }
  output_file=fopen(output_file_name,"w");
  fprintf(stdout,"\n\n");
  fprintf(stdout,"Length of file set to %d, target to %f.\n",
	  set_length,target_value);
  while(running_length<set_length){
    fscanf(trim_file,"%d\t%d\t%x\n",
	   &in_slot,&in_chan,&in_trim);
    fscanf(ped_file,"%d %f %d\n",
	   &in_ped_num,&in_ped,&in_thresh);
    /*    out_trim=((int)((target_value-in_ped)/ADC_PEDESTAL_SLOPE))+in_trim;*/
    out_trim=generate_new_trim(target_value,in_ped,in_trim);
    if(out_trim>0xff)
      out_trim=0xff;
    if(out_trim<0)
      out_trim=0;
    fprintf(output_file,"%d\t%d\t%x\n",in_slot,in_chan,out_trim);
    fprintf(stdout,"channel %d, ped=%6.3f, in_trim=%x, out_trim=%x\n",
	    running_length,in_ped,in_trim,out_trim);
    running_length++;
  }
  fprintf(stdout,"Finished generating new pedestal trim config file.\n");
  fclose(output_file);
  fclose(trim_file);
  fclose(ped_file);
}

int generate_new_trim(int target_value,
		      float current_pedestal,
		      int current_trim){
  float desired_delta_pedestal;
  float delta_trim_f;
  int delta_trim;
  int new_trim;
  
  desired_delta_pedestal=((float)target_value)-current_pedestal;
  /* new minus old ; the units are in ADC counts.  To convert the number to
     how far the pedestal trim has to move, we divide by the slope, 
     which is (ADC counts)/(trim count) */
  delta_trim_f=desired_delta_pedestal/ADC_PEDESTAL_SLOPE;
  delta_trim_f*=TRIM_CHANGE_THROTTLE;
  delta_trim=(int)delta_trim_f;
  new_trim=current_trim+delta_trim;
  return new_trim;
}
