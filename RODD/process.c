#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include "runs.h"
#include "errors.h"
#include "global.h"

#include <cernlib.h>

#include <makeHits.h>

/*  Hbook stuff */
#define HBOOK 800000
struct hbk{float ia[HBOOK];} pawc_ ;

static void display_lgd_raw(adc_values_t*);
static void display_rpd(adc_values_t*);
static void display_cpv(adc_values_t*);
static void display_upv(adc_values_t*);
static void display_bsd(adc_values_t*);
static void display_bgv(adc_values_t*);
static void display_tdc_rpd(tdc_values_t*);
static void display_tagger(tdc_values_t*);
static void display_cpv_tdcs(tdc_values_t*);
static void display_upv_tdcs(tdc_values_t*);
static void display_bgv_tdcs(tdc_values_t*);
static void display_bsd_tdcs(tdc_values_t*);
      
static adc_values_t *adc_vals, *rpd_adcs, *lgd_adcs, *cpv_adcs, *esum_adcs, *upv_adcs, *bsd_adcs, *bgv_adcs;
static tdc_values_t *tdc_vals, *rpd_tdcs, *tagger_tdcs, *cpv_tdcs, *upv_tdcs, *bgv_tdcs, *bsd_tdcs;
static scaler_values_t *scaler_vals;

/* 
 * At this point, all the tape data is read into memory and local storage;
 * user routines to do further data analysis and histogramming go here.
 */ 

int process(itape_header_t *event)
{
   int i;

   lgd_hits_t *lgdHits=NULL,*clusterHits=NULL;  
   lgd_clusters_t *clusters=NULL, *tmp=NULL;
 
   if(lgdHits=data_getGroup(event,GROUP_LGD_HITS,0)){
     hf1(3,lgdHits->nhits,1.0);
     for (i=0; i<lgdHits->nhits; i++){
       hf1(4,lgdHits->hit[i].energy,1.0);
     }

     if(clusters=data_getGroup(event,GROUP_LGD_CLUSTERS,0)){
       hf1(2,clusters->nClusters,1.0);
       for (i=0; i<clusters->nClusters; i++){
         hf1(5,clusters->cluster[i].nBlocks,1.0);
         hf1(6,clusters->cluster[i].width,1.0);
         hf1(7,clusters->cluster[i].energy,1.0);
         hf2(8,clusters->cluster[i].energy,clusters->cluster[i].width,1.0);
       }
     }
   }

   /*  Here we have simple functions to call user routines to histogram
    *      raw ADC and TDC data 
    *  This is where the user can insert code to do any higher-level
    *  analysis and/or histrogram the data in any way he/she chooses
    */

   /*  check that the trigger type is the correct one  */
   if (trig_flag==0 || trig_flag==event->trigger)  {

      hf1(1,event->trigger,1.0);

      if((lgd_adcs = data_getGroup(event,GROUP_LGD_ADCS,0)) != NULL){
        display_lgd_raw(lgd_adcs);
      }
      if((rpd_adcs = data_getGroup(event,GROUP_RPD_ADCS,0)) != NULL){
        display_rpd(rpd_adcs);  
      }
      if((cpv_adcs = data_getGroup(event,GROUP_CPV_ADCS,0)) != NULL){
        display_cpv(cpv_adcs);          
      }
      if((upv_adcs = data_getGroup(event,GROUP_UPV_ADCS,0)) != NULL){
        display_upv(upv_adcs);          
      }
      if((bsd_adcs = data_getGroup(event,GROUP_BSD_ADCS,0)) != NULL){
        display_bsd(bsd_adcs);          
      }
      if((bgv_adcs = data_getGroup(event,GROUP_BGV_ADCS,0)) != NULL){
        display_bgv(bgv_adcs);          
      }
      if((cpv_tdcs = data_getGroup(event,GROUP_CPV_TDCS,0)) != NULL){
        display_cpv_tdcs(cpv_tdcs);          
      }
      if((upv_tdcs = data_getGroup(event,GROUP_UPV_TDCS,0)) != NULL){
        display_upv_tdcs(upv_tdcs);          
      }
      if((bgv_tdcs = data_getGroup(event,GROUP_BGV_TDCS,0)) != NULL){
        display_bgv_tdcs(bgv_tdcs);          
      }
      if((bsd_tdcs = data_getGroup(event,GROUP_BSD_TDCS_LONG,0)) != NULL){
        display_bsd_tdcs(bsd_tdcs);          
      }
      if((rpd_tdcs = data_getGroup(event,GROUP_RPD_TDCS,0)) != NULL){
        display_tdc_rpd(rpd_tdcs); 
      }
      if((tagger_tdcs = data_getGroup(event,GROUP_TAGGER_TDCS,0)) != NULL){
        display_tagger(tagger_tdcs);
      }
   }
   return(0);
}

static void 
display_lgd_raw(adc_values_t *adc_Data) {
  int i, chan, val, icolumn, irow;
  int pedestal = 40;
  int dor, mam;

  for (i = 0; i < adc_Data->nadc; i++) {
    chan = adc_Data->adc[i].channel; val = adc_Data->adc[i].value; 
    dor = adc_Data->adc[i].dor;
    mam = adc_Data->adc[i].mam;


    /*
       Get the row and column number in the LGD from the software
          channel number ("chan")
    */
    irow=chan/28;
    icolumn = chan-(irow)*28;

    hf1(chan+10000,val,1.);
    hf1(chan+20000,val,1.);
    hf2(1000,icolumn,irow,val);
    hf1(11000+irow,icolumn,val);
    hf1(12000+icolumn,irow,val);

    if (dor) hf2(1002,icolumn,irow,1.);
    if (mam) hf2(1003,icolumn,irow,1.);

    if (val > pedestal) {
      hf1(999,chan,1.);
      hf2(1001,icolumn,irow,1.);
    }
  }
}


static void 
display_rpd(adc_values_t *adc_Data) {
  int i, chan, val;
  for (i = 0; i < adc_Data->nadc; i++) {
    chan = adc_Data->adc[i].channel; val = adc_Data->adc[i].value; 
    hf1(chan+2000,val,1.);
  }
}

static void 
display_cpv(adc_values_t *adc_Data) {
  int i, chan, val;
  for (i = 0; i < adc_Data->nadc; i++) {
    chan = adc_Data->adc[i].channel; val = adc_Data->adc[i].value; 
    hf1(chan+3000,val,1.);
 }
}
static void 
display_upv(adc_values_t *adc_Data) {
  int i, chan, val;
  for (i = 0; i < adc_Data->nadc; i++) {
    chan = adc_Data->adc[i].channel; val = adc_Data->adc[i].value; 
    hf1(chan+3500,val,1.);
  }
}
static void 
display_bsd(adc_values_t *adc_Data) {
  int i, chan, val;
  for (i = 0; i < adc_Data->nadc; i++) {
    chan = adc_Data->adc[i].channel; val = adc_Data->adc[i].value; 
    hf1(chan+6000,val,1.);
  }
}
static void 
display_bgv(adc_values_t *adc_Data) {
  int i, chan, val;
  for (i = 0; i < adc_Data->nadc; i++) {
    chan = adc_Data->adc[i].channel; val = adc_Data->adc[i].value; 
    hf1(chan+8000,val,1.);
  }
}

static void 
display_cpv_tdcs(tdc_values_t *tdc_Data) {
  int i, chan, val;
  for (i = 0; i < tdc_Data->ntdc; i++) {
    chan = tdc_Data->tdc[i].channel; val = tdc_Data->tdc[i].le; 
    hf1(chan+3700,val,1.);
 }
}

static void 
display_upv_tdcs(tdc_values_t *tdc_Data) {
  int i, chan, val;
  for (i = 0; i < tdc_Data->ntdc; i++) {
    chan = tdc_Data->tdc[i].channel; val = tdc_Data->tdc[i].le; 
    hf1(chan+7000,val,1.);
  }
}
static void 
display_bgv_tdcs(tdc_values_t *tdc_Data) {
  int i, chan, val;
  for (i = 0; i < tdc_Data->ntdc; i++) {
    chan = tdc_Data->tdc[i].channel; val = tdc_Data->tdc[i].le; 
    hf1(chan+8500,val,1.);
  }
}
static void 
display_bsd_tdcs(tdc_values_t *tdc_Data) {
  int i, chan, val;
  for (i = 0; i < tdc_Data->ntdc; i++) {
    chan = tdc_Data->tdc[i].channel; val = tdc_Data->tdc[i].le; 
    hf1(chan+6500,val,1.);
  }
}

static void 
display_tagger(tdc_values_t *tdc_Data) {
  int i, chan, val;
  for (i = 0; i < tdc_Data->ntdc; i++) {
    chan = tdc_Data->tdc[i].channel; val = tdc_Data->tdc[i].le; 
    hf1(chan+4000,val,1.);
    if (val>0 && val <1000) {
      if (chan < 19) {
      hf1(4100,chan,1.);
      }
      else {
      hf1(4200,chan-19,1.);    
      }
    }
  }
}

static void 
display_tdc_rpd(tdc_values_t *tdc_Data) {
  int i, chan, val;
  for (i = 0; i < tdc_Data->ntdc; i++) {
    chan = tdc_Data->tdc[i].channel; val = tdc_Data->tdc[i].le; 
    hf1(chan+4994,val,1.);
  }
}


 
void setupHbook(void)
     /*  Here's where we define all the histograms and scatterplots   */
{
    int i; 

    hlimit(HBOOK);

    hbook1(1,"trigger bit, all events",35,-0.5,34.5,0.0);
    hbook1(2,"number of clusters, physics events",12,-0.5,11.5,0.0);
    hbook1(3,"number of LGD hits, physics events",70,-0.5,69.5,0.0);
    hbook1(4,"energies of LGD hits, physics events",300,0.0,6.0,0.0);
    hbook1(5,"number blocks per cluster, physics events",20,0.0,20.0,0.0);
    hbook1(6,"width of clusters (cm), physics events",100,0.0,10.0,0.0);
    hbook1(7,"energies of clusters, physics events",300,0.0,6.0,0.0);
    hbook2(8,"width vs. energies of clusters, physics events",30,0.0,6.0,10,0.0,10.,0.0);

    /*
    for (i=10000; i<10784; i++) {
       hbook1(i," LGD channel ID-10000 ",200,0.0,4000.0,0.);
    }

    for (i=20000; i<20784; i++){
        hbook1(i," LGD channel ID-20000 ",100,-0.5,99.5,0.);
    }
    */

    hbook1(2000," RPD G0 ADC ",400,0.0,4000.0,0.);
    hbook1(2001," RPD G1 ADC ",400,0.0,4000.0,0.);
    hbook1(2002," RPD G2 ADC ",400,0.0,4000.0,0.);
    hbook1(2003," RPD G3 ADC ",400,0.0,4000.0,0.);
    hbook1(2004," RPD G4 ADC ",400,0.0,4000.0,0.);
    hbook1(2005," RPD G5 ADC ",400,0.0,4000.0,0.);
    hbook1(2006," RPD E0 ADC ",400,0.0,4000.0,0.);
    hbook1(2007," RPD E1 ADC ",400,0.0,4000.0,0.);
    hbook1(2008," RPD E2 ADC ",400,0.0,4000.0,0.);
    hbook1(2009," RPD E3 ADC ",400,0.0,4000.0,0.);
    hbook1(2010," RPD E4 ADC ",400,0.0,4000.0,0.);
    hbook1(2011," RPD E5 ADC ",400,0.0,4000.0,0.);
    hbook1(2012," RPD E6 ADC ",400,0.0,4000.0,0.);
    hbook1(2013," RPD E7 ADC ",400,0.0,4000.0,0.);
    hbook1(2014," RPD E8 ADC ",400,0.0,4000.0,0.);
    hbook1(2015," RPD E9 ADC ",400,0.0,4000.0,0.);
    hbook1(2016," RPD E10 ADC ",400,0.0,4000.0,0.);
    hbook1(2017," RPD E11 ADC ",400,0.0,4000.0,0.);

    hbook1(3000," CPV 0 ADC ",400,0.0,4000.0,0.);
    hbook1(3001," CPV 1 ADC ",400,0.0,4000.0,0.);
    hbook1(3002," CPV 2 ADC ",400,0.0,4000.0,0.);
    hbook1(3003," CPV 3 ADC ",400,0.0,4000.0,0.);
    hbook1(3004," CPV 4 ADC ",400,0.0,4000.0,0.);
    hbook1(3005," CPV 5 ADC ",400,0.0,4000.0,0.);
    hbook1(3006," CPV 6 ADC ",400,0.0,4000.0,0.);
    hbook1(3007," CPV 7 ADC ",400,0.0,4000.0,0.);
    hbook1(3008," CPV 8 ADC ",400,0.0,4000.0,0.);
    hbook1(3009," CPV 9 ADC ",400,0.0,4000.0,0.);
    hbook1(3010," CPV 10 ADC ",400,0.0,4000.0,0.);
    hbook1(3011," CPV 11 ADC ",400,0.0,4000.0,0.);
    hbook1(3012," CPV 12 ADC ",400,0.0,4000.0,0.);
    hbook1(3013," CPV 13 ADC ",400,0.0,4000.0,0.);
    hbook1(3014," CPV 14 ADC ",400,0.0,4000.0,0.);
    hbook1(3015," CPV 15 ADC ",400,0.0,4000.0,0.);
    hbook1(3016," CPV 16 ADC ",400,0.0,4000.0,0.);
    hbook1(3017," CPV 17 ADC ",400,0.0,4000.0,0.);
    hbook1(3018," CPV 18 ADC ",400,0.0,4000.0,0.);
    hbook1(3019," CPV 19 ADC ",400,0.0,4000.0,0.);
    hbook1(3020," CPV 20 ADC ",400,0.0,4000.0,0.);
    hbook1(3021," CPV 21 ADC ",400,0.0,4000.0,0.);
    hbook1(3022," CPV 22 ADC ",400,0.0,4000.0,0.);
    hbook1(3023," CPV 23 ADC ",400,0.0,4000.0,0.);
    hbook1(3024," CPV 24 ADC ",400,0.0,4000.0,0.);
    hbook1(3025," CPV 25 ADC ",400,0.0,4000.0,0.);
    hbook1(3026," CPV 26 ADC ",400,0.0,4000.0,0.);
    hbook1(3027," CPV 27 ADC ",400,0.0,4000.0,0.);
    hbook1(3028," CPV 28 ADC ",400,0.0,4000.0,0.);
    hbook1(3029," CPV 29 ADC ",400,0.0,4000.0,0.);

    hbook1(3500," UPV H0 ADC ",400,0.0,4000.0,0.);
    hbook1(3501," UPV H1 ADC ",400,0.0,4000.0,0.);
    hbook1(3502," UPV H2 ADC ",400,0.0,4000.0,0.);
    hbook1(3503," UPV H3 ADC ",400,0.0,4000.0,0.);
    hbook1(3504," UPV H4 ADC ",400,0.0,4000.0,0.);
    hbook1(3505," UPV H5 ADC ",400,0.0,4000.0,0.);
    hbook1(3506," UPV V0 ADC ",400,0.0,4000.0,0.);
    hbook1(3507," UPV V1 ADC ",400,0.0,4000.0,0.);

    hbook1(3700," CPV 0 TDC",300,0.0,300.0,0.);
    hbook1(3701," CPV 1 TDC",300,0.0,300.0,0.);
    hbook1(3702," CPV 2 TDC",300,0.0,300.0,0.);
    hbook1(3703," CPV 3 TDC",300,0.0,300.0,0.);
    hbook1(3704," CPV 4 TDC",300,0.0,300.0,0.);
    hbook1(3705," CPV 5 TDC",300,0.0,300.0,0.);
    hbook1(3706," CPV 6 TDC",300,0.0,300.0,0.);
    hbook1(3707," CPV 7 TDC",300,0.0,300.0,0.);
    hbook1(3708," CPV 8 TDC",300,0.0,300.0,0.);
    hbook1(3709," CPV 9 TDC",300,0.0,300.0,0.);
    hbook1(3710," CPV 10 TDC",300,0.0,300.0,0.);
    hbook1(3711," CPV 11 TDC",300,0.0,300.0,0.);
    hbook1(3712," CPV 12 TDC",300,0.0,300.0,0.);
    hbook1(3713," CPV 13 TDC",300,0.0,300.0,0.);
    hbook1(3714," CPV 14 TDC",300,0.0,300.0,0.);
    hbook1(3715," CPV 15 TDC",300,0.0,300.0,0.);
    hbook1(3716," CPV 16 TDC",300,0.0,300.0,0.);
    hbook1(3717," CPV 17 TDC",300,0.0,300.0,0.);
    hbook1(3718," CPV 18 TDC",300,0.0,300.0,0.);
    hbook1(3719," CPV 19 TDC",300,0.0,300.0,0.);
    hbook1(3720," CPV 20 TDC",300,0.0,300.0,0.);
    hbook1(3721," CPV 21 TDC",300,0.0,300.0,0.);
    hbook1(3722," CPV 22 TDC",300,0.0,300.0,0.);
    hbook1(3723," CPV 23 TDC",300,0.0,300.0,0.);
    hbook1(3724," CPV 24 TDC",300,0.0,300.0,0.);
    hbook1(3725," CPV 25 TDC",300,0.0,300.0,0.);
    hbook1(3726," CPV 26 TDC",300,0.0,300.0,0.);
    hbook1(3727," CPV 27 TDC",300,0.0,300.0,0.);
    hbook1(3728," CPV 28 TDC",300,0.0,300.0,0.);
    hbook1(3729," CPV 29 TDC",300,0.0,300.0,0.);

/*
    hbook1(3730," CPV OR TDC",500,0.0,1000.0,0.);
*/

    hbook1(4994," RPD G0 TDC ",100,0.0,1000.0,0.);
    hbook1(4995," RPD G1 TDC ",100,0.0,1000.0,0.);
    hbook1(4996," RPD G2 TDC ",100,0.0,1000.0,0.);
    hbook1(4997," RPD G3 TDC ",100,0.0,1000.0,0.);
    hbook1(4998," RPD G4 TDC ",100,0.0,1000.0,0.);
    hbook1(4999," RPD G5 TDC ",100,0.0,1000.0,0.);

    hbook1(5000," RPD E0 TDC ",100,0.0,1000.0,0.);
    hbook1(5001," RPD E1 TDC ",100,0.0,1000.0,0.);
    hbook1(5002," RPD E2 TDC ",100,0.0,1000.0,0.);
    hbook1(5003," RPD E3 TDC ",100,0.0,1000.0,0.);
    hbook1(5004," RPD E4 TDC ",100,0.0,1000.0,0.);
    hbook1(5005," RPD E5 TDC ",100,0.0,1000.0,0.);
    hbook1(5006," RPD E6 TDC ",100,0.0,1000.0,0.);
    hbook1(5007," RPD E7 TDC ",100,0.0,1000.0,0.);
    hbook1(5008," RPD E8 TDC ",100,0.0,1000.0,0.);
    hbook1(5009," RPD E9 TDC ",100,0.0,1000.0,0.);
    hbook1(5010," RPD E10 TDC ",100,0.0,1000.0,0.);
    hbook1(5011," RPD E11 TDC ",100,0.0,1000.0,0.);

    hbook1(5012," RPD OR TDC ",100,0.0,1000.0,0.);

 
    hbook1(4000," Tagger TDC L1 ",500,0.0,500.0,0.);
    hbook1(4001," Tagger TDC L2 ",500,0.0,500.0,0.);
    hbook1(4002," Tagger TDC L3 ",500,0.0,500.0,0.);
    hbook1(4003," Tagger TDC L4 ",500,0.0,500.0,0.);
    hbook1(4004," Tagger TDC L5 ",500,0.0,500.0,0.);
    hbook1(4005," Tagger TDC L6 ",500,0.0,500.0,0.);
    hbook1(4006," Tagger TDC L7 ",500,0.0,500.0,0.);
    hbook1(4007," Tagger TDC L8 ",500,0.0,500.0,0.);
    hbook1(4008," Tagger TDC L9 ",500,0.0,500.0,0.);
    hbook1(4009," Tagger TDC L10 ",500,0.0,500.0,0.);
    hbook1(4010," Tagger TDC L11 ",500,0.0,500.0,0.);
    hbook1(4011," Tagger TDC L12 ",500,0.0,500.0,0.);
    hbook1(4012," Tagger TDC L13 ",500,0.0,500.0,0.);
    hbook1(4013," Tagger TDC L14 ",500,0.0,500.0,0.);
    hbook1(4014," Tagger TDC L15 ",500,0.0,500.0,0.);
    hbook1(4015," Tagger TDC L16 ",500,0.0,500.0,0.);
    hbook1(4016," Tagger TDC L17 ",500,0.0,500.0,0.);
    hbook1(4017," Tagger TDC L18 ",500,0.0,500.0,0.);
    hbook1(4018," Tagger TDC L19 ",500,0.0,500.0,0.);
    hbook1(4019," Tagger TDC R1 ",500,0.0,500.0,0.);
    hbook1(4020," Tagger TDC R2 ",500,0.0,500.0,0.);
    hbook1(4021," Tagger TDC R3 ",500,0.0,500.0,0.);
    hbook1(4022," Tagger TDC R4 ",500,0.0,500.0,0.);
    hbook1(4023," Tagger TDC R5 ",500,0.0,500.0,0.);
    hbook1(4024," Tagger TDC R6 ",500,0.0,500.0,0.);
    hbook1(4025," Tagger TDC R7 ",500,0.0,500.0,0.);
    hbook1(4026," Tagger TDC R8 ",500,0.0,500.0,0.);
    hbook1(4027," Tagger TDC R9 ",500,0.0,500.0,0.);
    hbook1(4028," Tagger TDC R10 ",500,0.0,500.0,0.);
    hbook1(4029," Tagger TDC R11 ",500,0.0,500.0,0.);
    hbook1(4030," Tagger TDC R12 ",500,0.0,500.0,0.);
    hbook1(4031," Tagger TDC R13 ",500,0.0,500.0,0.);
    hbook1(4032," Tagger TDC R14 ",500,0.0,500.0,0.);
    hbook1(4033," Tagger TDC R15 ",500,0.0,500.0,0.);
    hbook1(4034," Tagger TDC R16 ",500,0.0,500.0,0.);
    hbook1(4035," Tagger TDC R17 ",500,0.0,500.0,0.);
    hbook1(4036," Tagger TDC R18 ",500,0.0,500.0,0.);
    hbook1(4037," Tagger TDC R19 ",500,0.0,500.0,0.);

    hbook1(4038," Tagger OR TDC  ",500,0.0,500.0,0.);

    hbook1(4100," Tagger Left occupancy ",19,-0.5,18.5,0.);
    hbook1(4200," Tagger Right occupancy ",19,-0.5,18.5,0.);


    hbook1(6000," BSD R0 ADC",400,0.0,4000.0,0.);
    hbook1(6001," BSD R1 ADC",400,0.0,4000.0,0.);
    hbook1(6002," BSD R2 ADC",400,0.0,4000.0,0.);
    hbook1(6003," BSD R3 ADC",400,0.0,4000.0,0.);
    hbook1(6004," BSD R4 ADC",400,0.0,4000.0,0.);
    hbook1(6005," BSD R5 ADC",400,0.0,4000.0,0.);
    hbook1(6006," BSD R6 ADC",400,0.0,4000.0,0.);
    hbook1(6007," BSD R7 ADC",400,0.0,4000.0,0.);
    hbook1(6008," BSD R8 ADC",400,0.0,4000.0,0.);
    hbook1(6009," BSD R9 ADC",400,0.0,4000.0,0.);
    hbook1(6010," BSD R10 ADC",400,0.0,4000.0,0.);
    hbook1(6011," BSD R11 ADC",400,0.0,4000.0,0.);
    hbook1(6012," BSD L0 ADC",400,0.0,4000.0,0.);
    hbook1(6013," BSD L1 ADC",400,0.0,4000.0,0.);
    hbook1(6014," BSD L2 ADC",400,0.0,4000.0,0.);
    hbook1(6015," BSD L3 ADC",400,0.0,4000.0,0.);
    hbook1(6016," BSD L4 ADC",400,0.0,4000.0,0.);
    hbook1(6017," BSD L5 ADC",400,0.0,4000.0,0.);
    hbook1(6018," BSD L6 ADC",400,0.0,4000.0,0.);
    hbook1(6019," BSD L7 ADC",400,0.0,4000.0,0.);
    hbook1(6020," BSD L8 ADC",400,0.0,4000.0,0.);
    hbook1(6021," BSD L9 ADC",400,0.0,4000.0,0.);
    hbook1(6022," BSD L10 ADC",400,0.0,4000.0,0.);
    hbook1(6023," BSD L11 ADC",400,0.0,4000.0,0.);
    hbook1(6024," BSD S0 ADC",400,0.0,4000.0,0.);
    hbook1(6025," BSD S1 ADC",400,0.0,4000.0,0.);
    hbook1(6026," BSD S2 ADC",400,0.0,4000.0,0.);
    hbook1(6027," BSD S3 ADC",400,0.0,4000.0,0.);
    hbook1(6028," BSD S4 ADC",400,0.0,4000.0,0.);
    hbook1(6029," BSD S5 ADC",400,0.0,4000.0,0.);
    hbook1(6030," BSD S6 ADC",400,0.0,4000.0,0.);
    hbook1(6031," BSD S7 ADC",400,0.0,4000.0,0.);
    hbook1(6032," BSD S8 ADC",400,0.0,4000.0,0.);
    hbook1(6033," BSD S9 ADC",400,0.0,4000.0,0.);
    hbook1(6034," BSD S10 ADC",400,0.0,4000.0,0.);
    hbook1(6035," BSD S11 ADC",400,0.0,4000.0,0.);
    hbook1(6036," BSD S12 ADC",400,0.0,4000.0,0.);
    hbook1(6037," BSD S13 ADC",400,0.0,4000.0,0.);
    hbook1(6038," BSD S14 ADC",400,0.0,4000.0,0.);
    hbook1(6039," BSD S15 ADC",400,0.0,4000.0,0.);
    hbook1(6040," BSD S16 ADC",400,0.0,4000.0,0.);
    hbook1(6041," BSD S17 ADC",400,0.0,4000.0,0.);
    hbook1(6042," BSD S18 ADC",400,0.0,4000.0,0.);
    hbook1(6043," BSD S19 ADC",400,0.0,4000.0,0.);
    hbook1(6044," BSD S20 ADC",400,0.0,4000.0,0.);
    hbook1(6045," BSD S21 ADC",400,0.0,4000.0,0.);
    hbook1(6046," BSD S22 ADC",400,0.0,4000.0,0.);
    hbook1(6047," BSD S23 ADC",400,0.0,4000.0,0.);

    hbook1(6500," BSD R0 TDC",300,0.0,300.0,0.);
    hbook1(6501," BSD R1 TDC",300,0.0,300.0,0.);
    hbook1(6502," BSD R2 TDC",300,0.0,300.0,0.);
    hbook1(6503," BSD R3 TDC",300,0.0,300.0,0.);
    hbook1(6504," BSD R4 TDC",300,0.0,300.0,0.);
    hbook1(6505," BSD R5 TDC",300,0.0,300.0,0.);
    hbook1(6506," BSD R6 TDC",300,0.0,300.0,0.);
    hbook1(6507," BSD R7 TDC",300,0.0,300.0,0.);
    hbook1(6508," BSD R8 TDC",300,0.0,300.0,0.);
    hbook1(6509," BSD R9 TDC",300,0.0,300.0,0.);
    hbook1(6510," BSD R10 TDC",300,0.0,300.0,0.);
    hbook1(6511," BSD R11 TDC",300,0.0,300.0,0.);
    hbook1(6512," BSD L0 TDC",300,0.0,300.0,0.);
    hbook1(6513," BSD L1 TDC",300,0.0,300.0,0.);
    hbook1(6514," BSD L2 TDC",300,0.0,300.0,0.);
    hbook1(6515," BSD L3 TDC",300,0.0,300.0,0.);
    hbook1(6516," BSD L4 TDC",300,0.0,300.0,0.);
    hbook1(6517," BSD L5 TDC",300,0.0,300.0,0.);
    hbook1(6518," BSD L6 TDC",300,0.0,300.0,0.);
    hbook1(6519," BSD L7 TDC",300,0.0,300.0,0.);
    hbook1(6520," BSD L8 TDC",300,0.0,300.0,0.);
    hbook1(6521," BSD L9 TDC",300,0.0,300.0,0.);
    hbook1(6522," BSD L10 TDC",300,0.0,300.0,0.);
    hbook1(6523," BSD L11 TDC",300,0.0,300.0,0.);
    hbook1(6524," BSD S0 TDC",300,0.0,300.0,0.);
    hbook1(6525," BSD S1 TDC",300,0.0,300.0,0.);
    hbook1(6526," BSD S2 TDC",300,0.0,300.0,0.);
    hbook1(6527," BSD S3 TDC",300,0.0,300.0,0.);
    hbook1(6528," BSD S4 TDC",300,0.0,300.0,0.);
    hbook1(6529," BSD S5 TDC",300,0.0,300.0,0.);
    hbook1(6530," BSD S6 TDC",300,0.0,300.0,0.);
    hbook1(6531," BSD S7 TDC",300,0.0,300.0,0.);
    hbook1(6532," BSD S8 TDC",300,0.0,300.0,0.);
    hbook1(6533," BSD S9 TDC",300,0.0,300.0,0.);
    hbook1(6534," BSD S10 TDC",300,0.0,300.0,0.);
    hbook1(6535," BSD S11 TDC",300,0.0,300.0,0.);
    hbook1(6536," BSD S12 TDC",300,0.0,300.0,0.);
    hbook1(6537," BSD S13 TDC",300,0.0,300.0,0.);
    hbook1(6538," BSD S14 TDC",300,0.0,300.0,0.);
    hbook1(6539," BSD S15 TDC",300,0.0,300.0,0.);
    hbook1(6540," BSD S16 TDC",300,0.0,300.0,0.);
    hbook1(6541," BSD S17 TDC",300,0.0,300.0,0.);
    hbook1(6542," BSD S18 TDC",300,0.0,300.0,0.);
    hbook1(6543," BSD S19 TDC",300,0.0,300.0,0.);
    hbook1(6544," BSD S20 TDC",300,0.0,300.0,0.);
    hbook1(6545," BSD S21 TDC",300,0.0,300.0,0.);
    hbook1(6546," BSD S22 TDC",300,0.0,300.0,0.);
    hbook1(6547," BSD S23 TDC",300,0.0,300.0,0.);

    hbook1(6548," BSD Right OR TDC",300,0.0,300.0,0.);
    hbook1(6549," BSD Left OR  TDC",300,0.0,300.0,0.);
    hbook1(6550," BSD Straight OR TDC",300,0.0,300.0,0.);
    hbook1(6551," BSD AND TDC",300,0.0,300.0,0.);

    hbook1(8000," BGV D0 ADC",400,0.0,4000.0,0.);
    hbook1(8001," BGV D1 ADC",400,0.0,4000.0,0.);
    hbook1(8002," BGV D2 ADC",400,0.0,4000.0,0.);
    hbook1(8003," BGV D3 ADC",400,0.0,4000.0,0.);
    hbook1(8004," BGV D4 ADC",400,0.0,4000.0,0.);
    hbook1(8005," BGV D5 ADC",400,0.0,4000.0,0.);
    hbook1(8006," BGV D6 ADC",400,0.0,4000.0,0.);
    hbook1(8007," BGV D7 ADC",400,0.0,4000.0,0.);
    hbook1(8008," BGV D8 ADC",400,0.0,4000.0,0.);
    hbook1(8009," BGV D9 ADC",400,0.0,4000.0,0.);
    hbook1(8010," BGV D10 ADC",400,0.0,4000.0,0.);
    hbook1(8011," BGV D11 ADC",400,0.0,4000.0,0.);
    hbook1(8012," BGV D12 ADC",400,0.0,4000.0,0.);
    hbook1(8013," BGV D13 ADC",400,0.0,4000.0,0.);
    hbook1(8014," BGV D14 ADC",400,0.0,4000.0,0.);
    hbook1(8015," BGV D15 ADC",400,0.0,4000.0,0.);
    hbook1(8016," BGV D16 ADC",400,0.0,4000.0,0.);
    hbook1(8017," BGV D17 ADC",400,0.0,4000.0,0.);
    hbook1(8018," BGV D18 ADC",400,0.0,4000.0,0.);
    hbook1(8019," BGV D19 ADC",400,0.0,4000.0,0.);
    hbook1(8020," BGV D20 ADC",400,0.0,4000.0,0.);
    hbook1(8021," BGV D21 ADC",400,0.0,4000.0,0.);
    hbook1(8022," BGV D22 ADC",400,0.0,4000.0,0.);
    hbook1(8023," BGV D23 ADC",400,0.0,4000.0,0.);
    hbook1(8024," BGV U0 ADC",400,0.0,4000.0,0.);
    hbook1(8025," BGV U1 ADC",400,0.0,4000.0,0.);
    hbook1(8026," BGV U2 ADC",400,0.0,4000.0,0.);
    hbook1(8027," BGV U3 ADC",400,0.0,4000.0,0.);
    hbook1(8028," BGV U4 ADC",400,0.0,4000.0,0.);
    hbook1(8029," BGV U5 ADC",400,0.0,4000.0,0.);
    hbook1(8030," BGV U6 ADC",400,0.0,4000.0,0.);
    hbook1(8031," BGV U7 ADC",400,0.0,4000.0,0.);
    hbook1(8032," BGV U8 ADC",400,0.0,4000.0,0.);
    hbook1(8033," BGV U9 ADC",400,0.0,4000.0,0.);
    hbook1(8034," BGV U10 ADC",400,0.0,4000.0,0.);
    hbook1(8035," BGV U11 ADC",400,0.0,4000.0,0.);
    hbook1(8036," BGV U12 ADC",400,0.0,4000.0,0.);
    hbook1(8037," BGV U13 ADC",400,0.0,4000.0,0.);
    hbook1(8038," BGV U14 ADC",400,0.0,4000.0,0.);
    hbook1(8039," BGV U15 ADC",400,0.0,4000.0,0.);
    hbook1(8040," BGV U16 ADC",400,0.0,4000.0,0.);
    hbook1(8041," BGV U17 ADC",400,0.0,4000.0,0.);
    hbook1(8042," BGV U18 ADC",400,0.0,4000.0,0.);
    hbook1(8043," BGV U19 ADC",400,0.0,4000.0,0.);
    hbook1(8044," BGV U20 ADC",400,0.0,4000.0,0.);
    hbook1(8045," BGV U21 ADC",400,0.0,4000.0,0.);
    hbook1(8046," BGV U22 ADC",400,0.0,4000.0,0.);
    hbook1(8047," BGV U23 ADC",400,0.0,4000.0,0.);

    hbook1(8500," BGV D0 TDC",500,0.0,500.0,0.);
    hbook1(8501," BGV D1 TDC",500,0.0,500.0,0.);
    hbook1(8502," BGV D2 TDC",500,0.0,500.0,0.);
    hbook1(8503," BGV D3 TDC",500,0.0,500.0,0.);
    hbook1(8504," BGV D4 TDC",500,0.0,500.0,0.);
    hbook1(8505," BGV D5 TDC",500,0.0,500.0,0.);
    hbook1(8506," BGV D6 TDC",500,0.0,500.0,0.);
    hbook1(8507," BGV D7 TDC",500,0.0,500.0,0.);
    hbook1(8508," BGV D8 TDC",500,0.0,500.0,0.);
    hbook1(8509," BGV D9 TDC",500,0.0,500.0,0.);
    hbook1(8510," BGV D10 TDC",500,0.0,500.0,0.);
    hbook1(8511," BGV D11 TDC",500,0.0,500.0,0.);
    hbook1(8512," BGV D12 TDC",500,0.0,500.0,0.);
    hbook1(8513," BGV D13 TDC",500,0.0,500.0,0.);
    hbook1(8514," BGV D14 TDC",500,0.0,500.0,0.);
    hbook1(8515," BGV D15 TDC",500,0.0,500.0,0.);
    hbook1(8516," BGV D16 TDC",500,0.0,500.0,0.);
    hbook1(8517," BGV D17 TDC",500,0.0,500.0,0.);
    hbook1(8518," BGV D18 TDC",500,0.0,500.0,0.);
    hbook1(8519," BGV D19 TDC",500,0.0,500.0,0.);
    hbook1(8520," BGV D20 TDC",500,0.0,500.0,0.);
    hbook1(8521," BGV D21 TDC",500,0.0,500.0,0.);
    hbook1(8522," BGV D22 TDC",500,0.0,500.0,0.);
    hbook1(8523," BGV D23 TDC",500,0.0,500.0,0.);
    hbook1(8524," BGV U0 TDC",500,0.0,500.0,0.);
    hbook1(8525," BGV U1 TDC",500,0.0,500.0,0.);
    hbook1(8526," BGV U2 TDC",500,0.0,500.0,0.);
    hbook1(8527," BGV U3 TDC",500,0.0,500.0,0.);
    hbook1(8528," BGV U4 TDC",500,0.0,500.0,0.);
    hbook1(8529," BGV U5 TDC",500,0.0,500.0,0.);
    hbook1(8530," BGV U6 TDC",500,0.0,500.0,0.);
    hbook1(8531," BGV U7 TDC",500,0.0,500.0,0.);
    hbook1(8532," BGV U8 TDC",500,0.0,500.0,0.);
    hbook1(8533," BGV U9 TDC",500,0.0,500.0,0.);
    hbook1(8534," BGV U10 TDC",500,0.0,500.0,0.);
    hbook1(8535," BGV U11 TDC",500,0.0,500.0,0.);
    hbook1(8536," BGV U12 TDC",500,0.0,500.0,0.);
    hbook1(8537," BGV U13 TDC",500,0.0,500.0,0.);
    hbook1(8538," BGV U14 TDC",500,0.0,500.0,0.);
    hbook1(8539," BGV U15 TDC",500,0.0,500.0,0.);
    hbook1(8540," BGV U16 TDC",500,0.0,500.0,0.);
    hbook1(8541," BGV U17 TDC",500,0.0,500.0,0.);
    hbook1(8542," BGV U18 TDC",500,0.0,500.0,0.);
    hbook1(8543," BGV U19 TDC",500,0.0,500.0,0.);
    hbook1(8544," BGV U20 TDC",500,0.0,500.0,0.);
    hbook1(8545," BGV U21 TDC",500,0.0,500.0,0.);
    hbook1(8546," BGV U22 TDC",500,0.0,500.0,0.);
    hbook1(8547," BGV U23 TDC",500,0.0,500.0,0.);


    hbook1(7000," UPV H0 TDC",1000,0.0,1000.0,0.);
    hbook1(7001," UPV H1 TDC",1000,0.0,1000.0,0.);
    hbook1(7002," UPV H2 TDC",1000,0.0,1000.0,0.);
    hbook1(7003," UPV H3 TDC",1000,0.0,1000.0,0.);
    hbook1(7004," UPV H4 TDC",1000,0.0,1000.0,0.);
    hbook1(7005," UPV H5 TDC",1000,0.0,1000.0,0.);
    hbook1(7006," UPV V0 TDC",1000,0.0,1000.0,0.);
    hbook1(7007," UPV V1 TDC",1000,0.0,1000.0,0.);

    hbook1(7008," UPV OR TDC",1000,0.0,1000.0,0.);

    /*
        scatterplot row and column, updated by actual ADC value 
    */
    hbook2(1000," Row vs. Column",28,-0.5,27.5,28,-0.5,27.5,0.);

    /*
        scatterplot row and column, number times hit above pedestal
    */
    hbook2(1001," Number hits above pedestal",28,-0.5,27.5,28,-0.5,27.5,0.);
    hbook2(1002," Number hits above DOR",28,-0.5,27.5,28,-0.5,27.5,0.);
    hbook2(1003," Number hits sent to MAM",28,-0.5,27.5,28,-0.5,27.5,0.);
     
    for (i=0; i<28; i++) {
      hbook1(11000+i," LGD Row+1100 = ID",28,-0.5,27.5,0.);
      hbook1(12000+i," LGD Column+1200 = ID",28,-0.5,27.5,0.);
    }
 
    hbook1(999," LGD channels above pedestal ",800,-0.5,799.5,0.);
}


void closeHbook()
{
    hrput(0,histo_file_name,"N");
}

