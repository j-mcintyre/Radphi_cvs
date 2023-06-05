#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calibration.h>
#include <lgdGeom.h>

#define N_LGD 784

static const char rcsid[]="$Id: ccToMap.c,v 1.1 1999/07/25 03:34:54 radphi Exp $";

void usage(void){
  fprintf(stdout,"usage of ccToMap:\n");
  fprintf(stdout,"ccToMap <ccfile> <run number>\n");
  fprintf(stdout,"<ccfile> is the cc.dat file from lgdCalibrate\n");
  fprintf(stdout,"<run number> is the run number where the ccs will take effect.\n\n");
}

void main(int nargs, char *args[])
{
  float cc[N_LGD];
  char in_char;
  int in_channel;
  float in_cc;
  int run_number=0;
  FILE *fp=NULL;
  int N_lines=0;
  int N_read;
  int i;
  for(i=0;i<N_LGD;i++){
    cc[i]=0.0;
  }
  if(nargs!=3){
    usage();
    exit(1);
  }
  fp = fopen(args[1],"r");
  if(!fp){
    fprintf(stderr,"Had a problem opening >%s< for input\n",args[1]);      
    exit(0);
  }
  run_number=atoi(args[2]);
  fprintf(stdout,"You would like to load LGD calibration constants from the\n");
  fprintf(stdout,"file named >%s< for run %d.  Is this correct? ",
	  args[1],run_number);
  fflush(stdout);
  in_char=getchar();
  if(in_char!='y'&&in_char!='Y'){
    fprintf(stdout,"Ok, not saving anything.\n");
    exit(0);
  }
  fprintf(stdout,"Now reading file >%s<.\n",args[1]);
  N_lines=0;
  while(!feof(fp)){
    fscanf(fp,"%f %d\n",&in_cc,&in_channel);
    if(in_channel<0||in_channel>N_LGD){
      fprintf(stderr,"There was a channel number out of bounds!\n");
      fprintf(stderr,"The line was:  %f %d\n",in_cc,in_channel);
      exit(1);
    }
    cc[in_channel]=in_cc;
    N_lines++;
  }
  fprintf(stdout,"Finished reading file, found %d lines.\n",N_lines);
  fprintf(stdout,"So you want to save the results of reading file\n");
  fprintf(stdout,">%s< as the cc result for run %d, correct? ",
	  args[1],run_number);
  fflush(stdout);
  fflush(stdin);
  in_char=getchar();
  if(in_char!='y'&&in_char!='Y'){
    fprintf(stdout,"Aborting save.\n");
    exit(0);
  }
  fprintf(stdout,"Ok, finally saving.\n");
  if(ccSave(CC_LGD,run_number,N_lines,cc)){
    fprintf(stderr,"Failed to save pedestals!\n");
    exit(1);
  }
  fprintf(stdout,"cc file uploaded to the map sucessfully.\n");
  exit(0);
}



