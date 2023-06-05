#include <stdio.h>
#include <stdlib.h>
#include <clib.h>
#include <math.h>

#define N_LGD_CHANNELS 784
#define LGD_TUBE_GAIN_EXP 9.362
/* this exponent stripped out of Scott's code. */
#define MAX_V 1950
#define MIN_V 1500
#define VOLTAGE_DUMMY_VALUE 0
#define CC_DUMMY_VALUE -.5

uint find_new_voltage(int V_old, float new_gain, float old_gain,
		      float exponent, int V_max, int V_min);
int read_cc_file(char *filename, float *cc, int *lines);
int read_voltage_file(char *filename, int *voltages, int *lines);
int user_didnt_type_yes(void);
void initialize_int_list(int *list, int value, int N);
void initialize_float_list(float *list, float value, int N);

uint find_new_voltage(int V_old, float new_gain, float old_gain,
		      float exponent, int V_max, int V_min){
  int newV;
  float log_of_q;
  float quotient;
  newV=((float)V_old)*exp(log(old_gain/new_gain)/exponent);
  if(newV>V_max)
    newV=V_max;
  if(newV<V_min)
    newV=V_min;
  return newV;
}

int read_cc_file(char *filename, float *cc, int *lines){
  FILE *infile=NULL;
  int in_channel;
  float in_cc;
  *lines=0;
  infile=fopen(filename,"r");
  if(infile==NULL)
    return(1);
  while(!feof(infile)){
    fscanf(infile,"%f %d\n",&in_cc,&in_channel);
    if(in_channel<0||in_channel>=N_LGD_CHANNELS)
      return(1);
    cc[in_channel]=in_cc;
    *lines++;
  }
  fclose(infile);
  return(0);
}

int read_voltage_file(char *filename, int *voltages, int *lines){
  FILE *infile=NULL;
  int in_channel,in_voltage;
  *lines=0;
  infile=fopen(filename,"r");
  if(infile==NULL)
    return(1);
  while(!feof(infile)){
    fscanf(infile,"%d %d\n",&in_channel,&in_voltage);
    if(in_channel<0||in_channel>=N_LGD_CHANNELS)
      return(1);
    voltages[in_channel]=in_voltage;
    *lines++;
  }
  fclose(infile);
  return(0);
}

void usage(void){
  fprintf(stdout,"setNewVoltages <voltage file> <cc file> <target cc>\n\n");
  return;
}

void initialize_int_list(int *list, int value, int N){
  int i;
  for(i=0;i<N;i++){
    list[i]=value;
  }
  return;
}

void initialize_float_list(float *list, float value, int N){
  int i;
  for(i=0;i<N;i++){
    list[i]=value;
  }
  return;
}

int user_didnt_type_yes(void){
  char in_char;
  fflush(stdout);
  fflush(stdin);
  in_char=getchar();
  if(in_char!='y'&&in_char!='Y'){
    return 1;
  }
  return 0;
}

void main(int nargs, char *args[]){
  FILE *outfile=NULL;
  char outfile_name[50]="new_voltages.dat";
  float cc[N_LGD_CHANNELS];
  float target_cc;
  int voltage[N_LGD_CHANNELS];
  int new_voltage;
  int i;
  int total;
  int lines;
  if(nargs!=4){
    usage();
    exit(1);
  }
  fprintf(stdout,"Voltage file=%s\n",args[1]);
  fprintf(stdout,"calibration constant file: %s\n",args[2]);
  target_cc=atof(args[3]);
  fprintf(stdout,"target calibration constant=%f\n",target_cc);
  fprintf(stdout,"new voltages will be written to: %s\n",outfile_name);
  fprintf(stdout,"Is this correct?");
  if(user_didnt_type_yes()){
    fprintf(stdout,"Exiting.\n");
    exit(0);
  }
  if((outfile=fopen(outfile_name,"w"))==NULL){
    fprintf(stderr,"A problem opening %s, exiting.\n",outfile_name);
  }
  initialize_int_list(voltage,VOLTAGE_DUMMY_VALUE,N_LGD_CHANNELS);
  initialize_float_list(cc,CC_DUMMY_VALUE,N_LGD_CHANNELS);
  if(read_voltage_file(args[1],voltage,&lines)){
    fprintf(stdout,"There was a problem loading %s, exiting.\n",args[1]);
    exit(1);
  }
  if(read_cc_file(args[2],cc,&lines)){
    fprintf(stdout,"There was a problem loading %s, exiting.\n",args[2]);
    exit(1);
  }
  total=0;
  for(i=0;i<N_LGD_CHANNELS;i++){
    if(cc[i]!=CC_DUMMY_VALUE && voltage[i]!=VOLTAGE_DUMMY_VALUE){
      new_voltage=(int)find_new_voltage((uint)voltage[i],target_cc,cc[i],
				   LGD_TUBE_GAIN_EXP,MAX_V,MIN_V);
      fprintf(outfile,"%d %d\n",i,new_voltage);
      total++;
    }
  }
  fprintf(stdout,"Generated %d new voltages.\n",total);
  fclose(outfile);
}
