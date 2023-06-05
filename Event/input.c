/* The event display here represented was designed, coded & documented by   */
/* Troels Petersen under the supervision of Elton Smith, with help from     */
/*   Dave  Heddle and James Muguira.					    */
/* Comments can be found around more complicated code leaving the basics    */
/* uncommented. Upon problems and/or questions, please contact Troels       */
/* Petersen at: <petersen@cebaf.gov> or <petersen@fys.ku.dk>		    */

#include "Event.h"

#include <disIO.h>
#include <disData.h>
#include <dataIO.h>
#include <itypes.h>
#include <iitypes.h>
#include <lgdCluster.h>
#include <eventType.h>
#include <triggerType.h>

#define BUFSIZE 100000

#define SQR(x) ((x)*(x))

int     fd;
uint32  *group_list;
itape_header_t *tape; 
int     mam_thresh;
char    *dispatcher=NULL;
char	command[1024];
Boolean fileopen=False;

/* - routine to select beam energy to get appropriate photon energy - */
/* - Lisa Kaufman - */

void BeamEnergySelect(Hv_Event hvevent)
{
  Widget dialog=NULL, amount;
  Widget rowcol, rc, menu;
  char   text[40];

  /* - creates a dialog box - */
  if(!dialog) {
    Hv_VaCreateDialog(&dialog, Hv_TITLE, "Select Beam Energy", NULL);
    rowcol = Hv_DialogColumn(dialog, 6);
    rc = Hv_DialogRow(rowcol, 6);
    amount = Hv_SimpleTextEdit(rc, "Beam Energy: (5.5=default)", NULL, 4);
    menu = Hv_StandardActionButtons(rowcol, 60, Hv_OKBUTTON + Hv_CANCELBUTTON);
  } /* end if */
  Hv_SetFloatText(amount, 5.5, 1);       /* sets default at 5.5 GeV */
  if(Hv_DoDialog(dialog, NULL) == Hv_OK){
    beamEnergy = Hv_GetFloatText(amount);  /* gets the value entered */
  } /* end if */
} /* end BeamEnergySelect */
  
/* - changes the color of an item (window comes up) - */

void EditColor(Hv_Event hvevent)
{
  Hv_Item	Item = hvevent->item;
  Hv_View	View = hvevent->view;
  short		color = Hv_GetItemColor(Item);

  Hv_ModifyColor(NULL, &color, "Change fill color", False);
  Hv_SetItemColor(Item, color);
  Hv_UpdateItem(Item);
}


/* - routine to go back one event - */

void ReadBackward(Hv_Event hvevent)
{
  int	i, j;
  
  if (fileopen) {
    if (event > 1) {
      --event;
      if (fd) close(fd); 
      fd = evt_open(filename, "r"); 
      for (i=0; i<event; i++) {
        evt_data_read(fd, tape, BUFSIZE);
      }
      MakeNew();
    }
    else {
      Hv_Information("No previous events!");
    }
  }
  else {
    Hv_Information("No data file is open!");
  }
}
  
/* - routine to read the next event - */

void ReadForward(Hv_Event hvevent)
{
  int		icount, j, ret, nread;
  int           running=1; 
  
  if(onlineflag==1) {

    while(running) {
      disIO_command("REQUEST_DISPLAY");
      ret = disIO_readCOOKED(tape,BUFSIZE,&nread,1);
      switch(ret){
      case DISIO_OK:
	/*         fprintf(stderr," calling GetData\n");  */

         GetData(tape,group_list);
         if  ( (trigflag!=trigbit)&(trigflag!=0) ){
           if (trigflag==100 & mam_energy >= mam_thresh){
	     MakeNew();
              return;
           }
         }
         if (trigflag==0){
	   MakeNew();
	   return;
         }
         else if (trigflag==trigbit){
	   MakeNew();
           return;
         }
         break;
      case DISIO_EOF:
         running=0;
         break;
      case DISIO_NODATA:
         break; 
      }
    }
  }

  if (fileopen) {
      if ((evt_data_read(fd, tape, BUFSIZE)) == DATAIO_OK) {
         GetData(tape,group_list);

         /* if what we just found was NOT a data event, we better
            try looking again */
         while (tape->eventType != EV_DATA) {
             if ((evt_data_read(fd, tape, BUFSIZE)) == DATAIO_OK) {
               GetData(tape,group_list); 
            } /* inner if*/
         } /* while */

         icount=0;
	 /*         fprintf(stderr,"trigflag=%d, trigbit=%d\n",trigflag,trigbit); */

         while  ( (trigflag!=trigbit)&(trigflag!=0) ){
           if (trigflag==100 & mam_energy >= mam_thresh){
              MakeNew();
              return;
           }
 	   if (icount <= BUFSIZE){
              icount++;
              if ((evt_data_read(fd, tape, BUFSIZE)) == DATAIO_OK) {  
                  GetData(tape,group_list);
		    /*   	          fprintf(stderr,"trigbit = %d; SKIPPED \n",trigbit); */
	      } else  {
               Hv_Information("At End-Of-File!");
	      }
           }
	   /* dsa - should probably pop up an info window here!! */
           else break;
 	 }
	 MakeNew();
    } else Hv_Information("At End-Of-File!");
  } else Hv_Information("No data file is open!");
}

/* - routine to select certain trigger types - */

void TrigSelect(Hv_Event hvevent)
{
  Widget dialog=NULL, amount;
  Widget rowcol, rc, menu;      
  char	 text[40];

  /* - creates a dialog box - */
  if (!dialog) {
    Hv_VaCreateDialog(&dialog, Hv_TITLE, "Select Trigger", NULL);
    rowcol = Hv_DialogColumn(dialog, 6);
    rc = Hv_DialogRow(rowcol, 6);
    amount = Hv_SimpleTextEdit(rc, "Trigger Type: (0=any) ", NULL, 4);
    menu = Hv_StandardActionButtons(rowcol, 60, Hv_OKBUTTON + Hv_CANCELBUTTON);
  }
  Hv_SetIntText(amount, trigval);		/* sets default to 0 */
  if (Hv_DoDialog(dialog, NULL) == Hv_OK)
    trigval = Hv_GetIntText(amount);	/* gets the value entered */
  switch(trigval){
  case 1: 
    if (trigflag==100) break;   /* if we have mam select, leave it on */
    trigflag=2; break;
  case 2: trigflag=4; break;
  case 3: trigflag=8; break;
  case 4: trigflag=16; break;
  case 5: trigflag=32; break; 
  default: trigflag=0;
  }
  ReadForward(hvevent);
}

/* - routine to select certain minimum MAM energy values  - */

void MAMSelect(Hv_Event	hvevent)
{
  Widget dialog=NULL, amount;
  Widget rowcol, rc, menu;      
  char	 text[40];

  /* - creates a dialog box - */
  if (!dialog) {
    Hv_VaCreateDialog(&dialog, Hv_TITLE, "Select mininum MAM Energy", NULL);
    rowcol = Hv_DialogColumn(dialog, 6);
    rc = Hv_DialogRow(rowcol, 6);
    amount = Hv_SimpleTextEdit(rc, "MAM Energy  ", NULL, 4);
    menu = Hv_StandardActionButtons(rowcol, 60, Hv_OKBUTTON + Hv_CANCELBUTTON);
  }
  Hv_SetIntText(amount, mam_thresh);		/* sets default to 0 */
  if (Hv_DoDialog(dialog, NULL) == Hv_OK){
    mam_thresh  = Hv_GetIntText(amount);	/* gets the value entered */
    trigflag=100;
  }
  ReadForward(hvevent);
}

/* - routine to jump to a certain event number - */

void ReadGoto(Hv_Event hvevent)
{
  int		i, j;
  Widget dialog=NULL, amount;
  Widget rowcol, rc, menu;      
  int	 value;
  char	 text[40];

  /* - creates a dialog box - */
  if (!dialog) {
    Hv_VaCreateDialog(&dialog, Hv_TITLE, "Go to event#", NULL);

    rowcol = Hv_DialogColumn(dialog, 6);
    rc = Hv_DialogRow(rowcol, 6);
    amount = Hv_SimpleTextEdit(rc, "Event number: ", NULL, 4);
    menu = Hv_StandardActionButtons(rowcol, 60, Hv_OKBUTTON + Hv_CANCELBUTTON);
  }
  Hv_SetIntText(amount, 1);		/* sets default to 1 */
  if (Hv_DoDialog(dialog, NULL) == Hv_OK)
    value = Hv_GetIntText(amount);	/* gets the value entered */

  /*  fprintf(stderr," onlineflag = %d\n",onlineflag);   */

  if (fileopen && (onlineflag==0)) {
    if (value < event) {
      if (fd) close(fd);
      fd = evt_open(filename, "r"); 
      event = 0;
    }
    while (evt_data_read(fd, tape, BUFSIZE) == DATAIO_OK) {
      GetData(tape,group_list);
      if (event >= value) {
        break;
      }
    }
    if (event == value) {
      MakeNew();
    }
    else {
      sprintf(text, "Search failed at event %d!", event);
      Hv_Information(text);
      MakeNew();
    }
  } else Hv_Information("No data file is open!");
}
  
/* - routine dealing with new data from standard input, not a file - */

void NewInput(Hv_Event hvevent)
{
  float angs, ange;
  int   xs,ys,xe,ye;
  int	i, j;
  char  text[100];
  int   channel; 
  int 	irow, icolumn; 
   
  if (fileopen) {
    if (fd) close(fd); 
  }
  event=0;

  tape = (itape_header_t *)malloc(BUFSIZE);
  xs=xc-(short)(210.0*cos(angs));
  ys=yc-(short)(210.0*sin(angs));
  xe=xc-(short)(210.0*cos(ange));
  ye=yc-(short)(210.0*sin(ange));
  Hv_SetLineWidth(2);
  Hv_DrawLine(xs, ys, xe, ye, Hv_black);
  group_list = (uint32 *)malloc(7*sizeof(uint32));

  fd = 0;

  fileopen=True;
  ReadForward(hvevent);
}

/* - routine dealing with new data from Dispatcher (online) - */

void Online(Hv_Event hvevent)
{
  int	index, i, j;
  char	text[100];
  int 	channel; 
  int	irow, icolumn; 
   
  if (fileopen) {
    close(fd); 
  }
  event=0;
 
  if(!dispatcher){
    dispatcher = getenv("RADPHI_DISPATCHER");
  }
  if(!dispatcher){
      sprintf(text, "No Dispatcher specified !\n");
      Hv_Information(text);
  }
    else { 
       if(disIO_connect(dispatcher,0))
          exit(1);

       sprintf(command,"Display:0x%x:0x%x",0xffff,0xffffffff);
       disIO_command(command);

       for(index=0;index<10;index++){
          disIO_command("REQUEST_DISPLAY");  
          }
    }

  tape = (itape_header_t *)malloc(BUFSIZE);
  group_list = (uint32 *)malloc(7*sizeof(uint32));
  fd = 0;
  onlineflag=1;
  
  fileopen=True;
  ReadForward(hvevent);

}

/* - routine dealing with new data - */

void NewFile(Hv_Event hvevent)
{
  int	i, j;
  char	text[100];
  int 	channel; 
  int 	irow, icolumn; 
   
  if (fileopen) {
    if (fd) close(fd); 
  }
  onlineflag=0;
  event=0;

  tape = (itape_header_t *)malloc(BUFSIZE);
  group_list = (uint32 *)malloc(7*sizeof(uint32));

  filename = Hv_FileSelect("Holy cow - real data!", "*.itape", NULL);
  fname = strrchr(filename,'/');
  ++fname;

  fd = evt_open(filename, "r"); 

  if (fd == 0)
    printf("Couldn't open the file (%s)!\n",filename);
  else {
    fileopen=True;
    ReadForward(hvevent);
  }
}

/* - routine dealing with itape data - */

void GetData(itape_header_t *tape, uint32 *list)
{
  int	i, j, k;
  char	text[100];

  adc_values_t *adcs;
  tdc_values_t *tdcs;
  reg_values_t *regs;
  scaler_values_t *scalers;
  int irow, icolumn;
  vector4_t pGamma[8];
  vector3_t vertex={0.0,0.0,0.0};

  lgd_hits_t *lgdHits=NULL;
  bsd_hits_t *bsdHits=NULL;
  bgv_hits_t *bgvHits=NULL;
  cpv_hits_t *cpvHits=NULL;
  upv_hits_t *upvHits=NULL;
  time_list_t *timelist=NULL;
  tagger_hits_t *taggerHits=NULL;
  photons_t *photons=NULL;
  mc_event_t *mc_event=NULL;
  mc_vertex_t *mc_vertex[10];
  reg_values_t *mam_regs = NULL;

  /* all these arrays and structs need to be initialized */
 
  mam_energy = 0;
  mam_mass = 0;
  scalers = NULL;
  
  for (i = 0; i < 6; i++){
    rpd_g[i] = 0;
  }
  for (i = 0; i < 12; i++){
    rpd_e[i] = 0;
  }
  for (i = 0; i < 28; i++){
    for (j = 0; j < 28; j++){
      lg[i][j] = 0;
    }
  }
  for (i = 0; i < 30; i++){
    cpv[i] = 0;
  }
  for (i = 0; i < 12; i++){
    bsd_r[i] = 0;
    bsd_l[i] = 0;
  }
  for (i = 0; i < 24; i++){
    bsd_s[i] = 0;
  }
  for (i = 0; i < 24; i++){
    bgv_ds[i] = 0;
    bgv_us[i] = 0;
  }
  for (i = 0; i < 8; i++){
    upv[i] = 0;
  }
  for (i = 0; i < 39; i++){
    tagger_times[i].multicount = 0;
  }
  for (i = 0; i < 12; i++){
    bsd_r_times[i].multicount = 0;
    bsd_l_times[i].multicount = 0;
  }
  for (i = 0; i < 24; i++){
    bsd_s_times[i].multicount = 0;
  }
  for (i = 0; i < 24; i++){
    bgv_ds_times[i].multicount = 0;
    bgv_us_times[i].multicount = 0;
  }
  for (i = 0; i < 30; i++){
    cpv_times[i].multicount = 0;
  }
  for (i = 0; i < 9; i++){
    upv_times[i].multicount = 0;
  }

  nhit=0;

  if( run != tape->runNo) {
    run = tape->runNo;
    if(setup_makeHits(tape)){
      fprintf(stderr,"Problem with setup_makeHits; exiting");
      return;
    }
  }
  
  if(tape->ngroups==1){
     unpackEvent(tape,tape->length);
  }

  /* If this isn't a data event, get out of here! */
  if (tape->eventType != EV_DATA) return;

  event =  tape->eventNo;
  trigbit = tape->trigger;

  /* We've got some data, so let's get higher-level groups!
     (if it is a physics event); */

  if ((trigbit == TRIG_DATA) || (trigbit == TRIG_MC)) {
    if (makeHits(tape) > 2) {
       fprintf(stdout, "problem with makeHits event %d\n", event);
       return;
    }
  }

  calibrated = pedestals;

  timelist = data_getGroup(tape,GROUP_TIME_LIST,0);

  if (pedestals) {
    if (lgdHits = data_getGroup(tape, GROUP_LGD_HITS, 0)) {
      for (j = 0; j < lgdHits->nhits; j++) {
	int irow = lgdHits->hit[j].channel / 28 ; 
	int icolumn = lgdHits->hit[j].channel - (28*irow) ; 
	lg[icolumn][irow] = (lgdHits->hit[j].energy)*1000.; /* GeV -> MeV */
	/*
	 fprintf(stderr,"lg value %d\n", (int)lg[icolumn][irow]);
	 fprintf(stderr,"irow %d\n", irow);
	 fprintf(stderr,"icolumn %d\n", icolumn);
	 fprintf(stderr,"energy value  %f\n", lgdHits->hit[j].energy);
	*/
      }
    }
  }
  else {
    if (adcs = data_getGroup(tape, GROUP_LGD_ADCS, 0)) {
      for (j = 0; j < adcs->nadc; j++) {
        int channel = adcs->adc[j].channel;
        int irow = channel / 28 ;
        int icolumn = channel - (28*irow) ; 
        lg[icolumn][irow] = adcs->adc[j].value;
      }
    }
  }

 if (adcs = data_getGroup(tape, GROUP_RPD_ADCS, 0)) {
    for (j = 0; j < adcs->nadc; j++) {
      int channel = adcs->adc[j].channel;
      float value = adcs->adc[j].value;
      if (channel < 6) {
        rpd_g[channel] = value;
      } 
      else {
        rpd_e[channel-6] = value;
      }
      /* fprintf(stdout, "RPD_ADCS:\n");  */
    }
  }
  if (tdcs = data_getGroup(tape, GROUP_RPD_TDCS, 0)) {
    /* fprintf(stdout, "RPD_TDCS:\n");  */
  }

  if (pedestals) {
    if (bsdHits = data_getGroup(tape,GROUP_BSD_HITS,0)) {
      for (j = 0; j < bsdHits->nhits; j++) {
        int channel = bsdHits->hit[j].channel;
        float value = (bsdHits->hit[j].energy)*1e6; /* in KeV */
        /*
         fprintf(stderr,"BSD channel  %d\n",  bsdHits->hit[j].channel);
         fprintf(stderr,"energy value  %f\n", bsdHits->hit[j].energy);
         fprintf(stderr,"value  %d\n", (int)value);
        */
        if (channel < 12) {
          bsd_r[channel] = value;
        }
        else if (channel < 24) {
          bsd_l[channel-12] = value;
        }
        else if (channel < 48) {
          bsd_s[channel-24] = value;
        }
        for (k = 0; k < bsdHits->hit[j].times; k++) {
          if (channel < 12){ 
            bsd_r_times[channel].time[bsd_r_times[channel].multicount++] =
                                 timelist->element[bsdHits->hit[j].tindex +k].le;
          }
          else if (channel < 24) {
            bsd_l_times[channel-12].time[bsd_l_times[channel-12].multicount++] =
                                    timelist->element[bsdHits->hit[j].tindex +k].le;
          }
          else {
            bsd_s_times[channel-24].time[bsd_s_times[channel-24].multicount++] =
                                    timelist->element[bsdHits->hit[j].tindex +k].le;
          }
  	}
      }
    }
  }
  else {
    if (adcs = data_getGroup(tape, GROUP_BSD_ADCS, 0)) {
      for (j = 0; j < adcs->nadc; j++) {
        int channel = adcs->adc[j].channel;
        float value = adcs->adc[j].value;
        if (channel < 12) {
  	  bsd_r[channel] = value;
        }
        else if (channel < 24) {
          bsd_l[channel-12] = value;
        }
        else if (channel < 48) {
          bsd_s[channel-24] = value;
        }
        /* fprintf(stdout, "BSD_ADCS:\n"); */
      }
    }
    if (tdcs = data_getGroup(tape, GROUP_BSD_TDCS_LONG, 0)) {
      for (j = 0; j < tdcs->ntdc; j++) {
        int channel = tdcs->tdc[j].channel;
        float value = tdcs->tdc[j].le/2.;
        if (channel < 12) { 
  	  bsd_r_times[channel].time[bsd_r_times[channel].multicount++] = value;
        }
        else if (channel < 24) {
  	  bsd_l_times[channel-12].time[bsd_l_times[channel-12].multicount++] = value;
        }
        else if (channel < 48) {
          bsd_s_times[channel-24].time[bsd_s_times[channel-24].multicount++] = value;
        }
      }
    }
  }

  if (pedestals) {
    if (bgvHits = data_getGroup(tape,GROUP_BGV_HITS,0)) {
      for (j = 0; j < bgvHits->nhits; j++) {
        int channel = bgvHits->hit[j].channel;
	bgv_us[channel] = bgvHits->hit[j].energy[0]*1000.;
	bgv_ds[channel] = bgvHits->hit[j].energy[1]*1000.;
        for (k = 0; k < bgvHits->hit[j].times[0]; k++) {
          bgv_us_times[channel].time[bgv_us_times[channel].multicount++] =
                                timelist->element[bgvHits->hit[j].tindex[0] +k].le;
        }
        for (k = 0; k < bgvHits->hit[j].times[0]; k++) {
          bgv_ds_times[channel].time[bgv_ds_times[channel].multicount++] =
                                timelist->element[bgvHits->hit[j].tindex[1] +k].le;
        }
      }
    }
  }
  else {
    if (adcs = data_getGroup(tape, GROUP_BGV_ADCS, 0)) {
      for (j = 0; j < adcs->nadc; j++) {
        int channel = adcs->adc[j].channel;
        float value = adcs->adc[j].value;
        if (channel < 24) {
          bgv_ds[channel] = value;
        }
        else {
          bgv_us[channel-24] = value;
        }
        /* fprintf(stdout, "BGV_ADCS:\n"); */
      }
    }
    if (tdcs = data_getGroup(tape, GROUP_BGV_TDCS, 0)) {
      for (j = 0; j < tdcs->ntdc; j++) {
        int channel = tdcs->tdc[j].channel;
        float value = tdcs->tdc[j].le/2.;
        if (channel < 24){
          bgv_ds_times[channel].time[bgv_ds_times[channel].multicount++] = value;
        }
        else {
          bgv_us_times[channel-24].time[bgv_us_times[channel-24].multicount++] = value;
        }
      }
    }
  }

  if (pedestals) {
    if (cpvHits = data_getGroup(tape,GROUP_CPV_HITS,0)) {
      for (j = 0; j < cpvHits->nhits; j++) {
        int channel = cpvHits->hit[j].channel;
        cpv[channel] = (cpvHits->hit[j].energy)*1e6;
	/*
	 fprintf(stderr,"CPV channel  %d\n",  cpvHits->hit[j].channel);
	 fprintf(stderr,"energy value  %f\n", cpvHits->hit[j].energy);
	 fprintf(stderr,"value  %d\n", (int)cpv[channel]);
	*/
        for (k = 0; k < cpvHits->hit[j].times; k++) {
          cpv_times[channel].time[cpv_times[channel].multicount++] =
                             timelist->element[cpvHits->hit[j].tindex +k].le;
        }
      }
    }
  }
  else {
    if (adcs = data_getGroup(tape, GROUP_CPV_ADCS, 0)) {
      for (j = 0; j < adcs->nadc; j++) {
        int channel = adcs->adc[j].channel;
        cpv[channel] = adcs->adc[j].value;
      }
      /* fprintf(stdout, "CPV_ADCS:\n");  */
    }
    if (tdcs = data_getGroup(tape, GROUP_CPV_TDCS, 0)) {
      for (j = 0; j < tdcs->ntdc; j++) {
        int channel = tdcs->tdc[j].channel;
        float value = tdcs->tdc[j].le/2.;
        if (channel < 30){
  	  cpv_times[channel].time[cpv_times[channel].multicount++] = value;
        }
      }
    }
  }

  if (pedestals) {
    if (upvHits = data_getGroup(tape,GROUP_UPV_HITS,0)) {
      for (j = 0; j < upvHits->nhits; j++) {
        int channel = upvHits->hit[j].channel;
        upv[channel] = (upvHits->hit[j].energy)*1e6;
        /*
         fprintf(stderr,"UPV channel  %d\n",  upvHits->hit[j].channel);
         fprintf(stderr,"energy value  %f\n", upvHits->hit[j].energy);
         fprintf(stderr,"value  %d\n", (int)upv[channel]);
        */
        for (k = 0; k < upvHits->hit[j].times; k++) {
          upv_times[channel].time[upv_times[channel].multicount++] =
                             timelist->element[upvHits->hit[j].tindex +k].le;
        }
      }
    }
  }
  else {
    if (adcs = data_getGroup(tape, GROUP_UPV_ADCS, 0)) {
      for (j = 0; j < adcs->nadc; j++) {
        int channel = adcs->adc[j].channel;
        upv[channel] = adcs->adc[j].value;
      }
    }
    if (tdcs = data_getGroup(tape, GROUP_UPV_TDCS, 0)) {
      for (j = 0; j < tdcs->ntdc; j++) {
        int channel = tdcs->tdc[j].channel;
        float value = tdcs->tdc[j].le/2.;
        if (channel < 8){
          upv_times[channel].time[upv_times[channel].multicount++] = value;
        }
      }
    }
  }

  if (pedestals) {
    if (taggerHits = data_getGroup(tape,GROUP_TAGGER_HITS,0)) {
      for (j = 0; j < taggerHits->nhits; j++) {
        int channel = taggerHits->hit[j].channel;
        /*
         fprintf(stderr,"tagger channel  %d\n",  taggerHits->hit[j].channel);
         fprintf(stderr,"energy value  %f\n", taggerHits->hit[j].energy);
         fprintf(stderr,"time  %d\n", taggerHits->hit[j].times);
        */
        for (k = 0; k < taggerHits->hit[j].times[0]; k++) {
          tagger_times[channel].time[tagger_times[channel].multicount++] =
                   timelist->element[taggerHits->hit[j].tindex[0] +k].le;
        }
        for (k = 0; k < taggerHits->hit[j].times[1]; k++) {
          tagger_times[channel+19].time[tagger_times[channel+19].multicount++] =
                   timelist->element[taggerHits->hit[j].tindex[1] +k].le;
        }
      }
    }
  }
  else {
    if (tdcs = data_getGroup(tape, GROUP_TAGGER_TDCS, 0)) {
      for (j = 0; j < tdcs->ntdc; j++) {
        int channel = tdcs->tdc[j].channel; 
        int tdc = tdcs->tdc[j].le/2.;
        tagger_times[channel].time[tagger_times[channel].multicount++] = tdc;
        /*  fprintf(stdout," tagger channel %d, value = %d\n", channel,tagger[channel] ); */
      }
      /* fprintf(stdout, "TAGGER_TDCS:\n"); */ 
    }
  }

  if (scalers = data_getGroup(tape, GROUP_SCALERS, 0)) {
    /*    removed for now - we should build in real scaler handling   */
    /*         fprintf(stdout, "--------SCALER data found:---------\n");  */
    /*         for (j=0; j < scaler_vals->nscalers; j++) {                */
    /*           fprintf(stdout, "  scaler channel: %d   value:  %d\n",j,scaler_vals->scaler[j]); */
    /*         } */
  }

  if (mam_regs = data_getGroup(tape, GROUP_MAM_REGS, 0)) {
    mam_energy = mam_regs->reg[7].contents;
    mam_mass = mam_regs->reg[8].contents;
  }

  noMC = 0;
  if (mc_event =  data_getGroup(tape, GROUP_MC_EVENT, 0)) {
    int noVert = mc_event->nfinalVert;	
    mcEventPrint(mc_event);
    for (j=0; j < noVert; j++){
      float zimp = 120;
      int l;
      mc_vertex[j] = data_getGroup(tape,GROUP_MC_EVENT,j+1);
      mcVertexPrint(noVert, mc_vertex[j]);
      for (l=0; l < mc_vertex[j]->npart; l++) {
        pxMC[noMC] = mc_vertex[j]->part[l].momentum.space.x;
        pyMC[noMC] = mc_vertex[j]->part[l].momentum.space.y;
        pzMC[noMC] = mc_vertex[j]->part[l].momentum.space.z;
        eMC[noMC] = sqrt(pxMC[noMC]*pxMC[noMC]
                        +pyMC[noMC]*pyMC[noMC]
                        +pzMC[noMC]*pzMC[noMC]);
        xMC[noMC] = zimp*(pxMC[noMC]/eMC[noMC]);
        yMC[noMC] = zimp*(pyMC[noMC]/eMC[noMC]);
        zMC[noMC] = zimp;
        noMC++;
      } 
    }
    noMC -= mc_event->npart;  /* last particle is recoil nucleon */
  }

  if (photons = data_getGroup(tape,GROUP_PHOTONS,0)) {
    nhit = photons->nPhotons;
    for (j=0; j<nhit; j++){
      e[j] = photons->photon[j].energy;
      px[j] = photons->photon[j].momentum.x;
      py[j] = photons->photon[j].momentum.y;
      pz[j] = photons->photon[j].momentum.z;
      x[j] = photons->photon[j].position.x;
      y[j] = photons->photon[j].position.y;
      z[j] = photons->photon[j].position.z;
    }
  }
}

/* prints the current mc_event group */

int mcEventPrint (mc_event_t *mc_event)
{

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
  printf("   p:  %f\n",sqrt(SQR(mc_event->part[k].momentum.space.x)
                           +SQR(mc_event->part[k].momentum.space.y)
                           +SQR(mc_event->part[k].momentum.space.z)));
  printf(" Printing Final State\n");
  printf("  Number of vertices: %i\n",mc_event->nfinalVert);
 }
 return 1; 
}

/* prints the current mc_vertex group */

int mcVertexPrint (int iv, mc_vertex_t *mc_vertex)
{

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
   printf("     p:   %f\n",sqrt(SQR(mc_vertex->part[l].momentum.space.x)
                               +SQR(mc_vertex->part[l].momentum.space.y)
                               +SQR(mc_vertex->part[l].momentum.space.z)));
 } 
 printf("------------------End of MC_GROUP--------------------------\n");
 return 1; 
}
