/* The event display here represented was designed, coded & documented by   */
/* Troels Petersen under the supervision of Elton Smith, with help from     */
/*   Dave  Heddle and James Muguira.					    */
/* Comments can be found around more complicated code leaving the basics    */
/* uncommented. Upon problems and/or questions, please contact Troels       */
/* Petersen at: <petersen@cebaf.gov> or <petersen@fys.ku.dk>		    */

#include "Event.h"
#include <triggerType.h>
#include <detector_params.h>
#include <makeHits.h>

void MakeNew(void)
{
/* - routine dealing with new data                                  -*/
/* - the routine also generates raw data for display and this part - */
/* - might therefore be worth keeping for testing the display. -     */
/*                                                                   */
/*  most of this commented out or removed (DSA) for real data        */
/*   also, initialization of mark[] changed to 8 values from 10      */

  Hv_Rect	rect, area;
  int		i, j, k, hit=1, length;
  float		er=7000.0, temp, dist;	
  float		theta, phi;	
  
  for (i=0; i < 8; i++) mark[i]=0;
  
  /* - updates the display with the new event - */
  Hv_DrawItem(LeadGlassItem, NULL);
  Hv_DrawItem(TaggerItem, NULL);
  Hv_DrawItem(CPVItem, NULL);
  Hv_DrawItem(BSDItem, NULL);
  Hv_DrawItem(UPVItem, NULL);
  Hv_DrawItem(Box4, NULL);
  Hv_CopyRect(&area, Hv_GetItemArea(Box5));
  Hv_SetRect(&rect, area.left+3, area.top+3, area.width-6, area.height-6);
  Hv_FillRect(&rect, Hv_gray12);
}

void LeadGlass(Hv_Item Item, Hv_Region region)
{
  /* - the area of the item and a dummy area for drawing the LGB's color - */
  Hv_Rect	area, rect;
  Hv_Rect       drawArea;

  /* - the radius of the photon circle - */
  int		rad;

  int		xs, ys;
  int		i, j;
  char		text[40];
  Hv_String	*Hvtext;
  int          	icolor;
  float		factor;

  Hv_CopyRect(&area, Hv_GetItemArea(Item));     	/* gets the area */
  drawArea = area;
  drawArea.left -= 10;
  drawArea.top -= 10;
  drawArea.width += 10;
  drawArea.height += 10;
  Hv_FillRect(&drawArea, Hv_gray12);			/* clears the area */
  Hv_GetItemCenter(Item, &xc, &yc);
  xc += -20;						/* moves the center */
  yc += 0;						/* and thereby the item */
  Hv_SetLineWidth(3);
  Hv_FillCircle(xc, yc, 180, True, Hv_white, Hv_black);
  Hv_SetLineWidth(1);
  xp = xc - 14*nsize;					/* finds the upper left corner */
  yp = yc - 14*nsize;

  
  if (strlen(fname)!=0){
    sprintf(text, "Filename: %s   Run: %d", fname, run); 
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left, area.top, Hvtext, NULL);
  }

  sprintf(text, "Event number: %d", event);
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left, area.top+12, Hvtext, NULL);

  if (trigbit==TRIG_DATA){
    sprintf(text, "Physics Trigger");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+300, area.top, Hvtext, NULL);
    if (mam_display==1) {
       sprintf(text, "MAM Energy: %d",mam_energy);
       Hvtext = Hv_CreateString(text);		   
       Hv_StringDraw(area.left+200, area.top+470, Hvtext, NULL);

       sprintf(text, "MAM Mass squared: %d",mam_mass);
       Hvtext = Hv_CreateString(text);		   
       Hv_StringDraw(area.left+200, area.top+480, Hvtext, NULL);
       }
  }
  else if (trigbit==TRIG_LGDMON){
    sprintf(text, "LGD Monitor Trigger");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+300, area.top, Hvtext, NULL);
  }
  else if (trigbit==TRIG_BGVMON){
    sprintf(text, "BGV Monitor Trigger");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+300, area.top, Hvtext, NULL);
  }
  else if (trigbit==TRIG_PED){
    sprintf(text, "Pedestal Trigger");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+300, area.top, Hvtext, NULL);
  }
  else if (trigbit==TRIG_BT){
    sprintf(text, "Basetest/scaler Trigger");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+300, area.top, Hvtext, NULL);
  }
  else {
  sprintf(text, "Unknown Trigger");
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left+300, area.top, Hvtext, NULL);
  }

  sprintf(text, "(Looking Downstream)");
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left, area.top+450, Hvtext, NULL);
  Hv_DestroyString(Hvtext);

  if (pedestals == 1 && (!calibrated)) {
    sprintf(text, "Warning:  Uncalibrated Data!");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+30, area.top+24, Hvtext, NULL);
    Hv_DestroyString(Hvtext);
  }
  else if (pedestals == 0 && calibrated) {
    sprintf(text, "Warning: no Pedestal Subtraction!");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+30, area.top+30, Hvtext, NULL);
    Hv_DestroyString(Hvtext);
  }
  else if (pedestals == 1 && (!calibrated)){
    sprintf(text, "Warning: Uncalibrated Data!");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+30, area.top+30, Hvtext, NULL);
    Hv_DestroyString(Hvtext);
  }
  else if (trigbit == TRIG_DATA){
    sprintf(text, "Calibrated data");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+30, area.top+30, Hvtext, NULL);
    Hv_DestroyString(Hvtext);
  }
  if (monte_carlo){
    sprintf(text, "Monte Carlo data");
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+30, area.top+43, Hvtext, NULL);
    Hv_DestroyString(Hvtext);
  }

  /* - loop for coloring the LGBs - */

  factor=(Hv_red-Hv_white)/(float)LGD_RED_HOT;    /* color scaling factor */
  if (trigbit==TRIG_BT) {
    factor=(Hv_red-Hv_white)/50.;   /* compress scale for basetest events */ 
  }

  for (i=0; i < 28; i++) {
    for (j=0; j < 28; j++) {
      if (lg[i][j] != 0) {
        float energy = lg[i][j];
	xs = xp + nsize*(27-i);
	ys = yp + nsize*(27-j);
	Hv_SetRect(&rect, xs, ys, nsize, nsize);
	if (calibrated) {
          if (energy > lb_lgdE) {
	    icolor = energy*factor + Hv_white;
	    if (icolor > Hv_white) icolor= Hv_white;
	    if (icolor < Hv_red) icolor= Hv_red;
	    Hv_FillRect(&rect, icolor);   /* scale: white - red */
	  }
	  else {
            Hv_FillRect(&rect, Hv_white); /* pedestal or below  */  
          }
	}
	else if (lg[i][j] > LB_ADC) {
	  icolor= lg[i][j]*factor + Hv_white;
	  if (icolor > Hv_white) icolor=Hv_white;
	  if (icolor < Hv_red) icolor=Hv_red;
	  Hv_FillRect(&rect, icolor);   /* scale: white - red */
	}
      }
    }
  }

  /* - draw the lines separating the LGBs - */

  for (i=0; i<15; i++) {
    Hv_DrawLine(xp+nsize*i, yc-nsize*geo[i], xp+nsize*i, yc+nsize*geo[i], Hv_black);
    Hv_DrawLine(xp+nsize*(28-i), yc-nsize*geo[i], xp+nsize*(28-i), yc+nsize*geo[i], Hv_black);
    Hv_DrawLine(xc-nsize*geo[i], yp+nsize*i, xc+nsize*geo[i], yp+nsize*i, Hv_black);
    Hv_DrawLine(xc-nsize*geo[i], yp+nsize*(28-i), xc+nsize*geo[i], yp+nsize*(28-i), Hv_black);
  }

  /* - mark the dead blocks with a cross - */

  for (i=0; i<lgd_nChannels; i++) {
    int row=i/lgd_nRows;
    int col=i%lgd_nRows;
    int halfwidth=geo[(row<14)? row:27-row];
    if ((channel_status_lgd[i] != CHANNEL_STATUS_GOOD) &&
        (col > 14-halfwidth) && (col < 13+halfwidth)) {
      vector3_t center;
      float x[2],y[2];
      lgdLocalCoord(i,&center);
      x[0]=center.x-2.0;
      x[1]=center.x+2.0;
      y[0]=center.y-2.0;
      y[1]=center.y+2.0;
      Hv_DrawLine(xc-(x[0]*ppc), yc-(y[0]*ppc), xc-(x[1]*ppc), yc-(y[1]*ppc), Hv_black);
      Hv_DrawLine(xc-(x[0]*ppc), yc-(y[1]*ppc), xc-(x[1]*ppc), yc-(y[0]*ppc), Hv_black);
    }
  }

  /* - draw beamhole and target - */

  Hv_SetRect(&rect, xc-nsize+1, yc-nsize+1, 2*nsize-1, 2*nsize-1);
  Hv_FillRect(&rect, Hv_white);
  Hv_SetRect(&rect, xc-3, yc-3, 7, 7);
  Hv_FillRect(&rect, Hv_black);

  /* - draw circle (rad = energy) and marks for each photon - */

  Hv_SetLineWidth(2);
  for (i=0; i < nhit; i++) {
    xs = xc - (int)(x[i]*ppc);
    ys = yc - (int)(y[i]*ppc);
    rad = (int)(e[i]*12.0);
    Hv_FrameCircle(xs, ys, rad, Hv_black);
    
    sprintf(text, "%d", i+1);
    Hvtext = Hv_CreateString(text);		   
    Hvtext->font=Hv_timesBold17;   
    Hv_StringDraw(xs-rad*0.8-8, ys-rad*0.8-8, Hvtext, NULL);
    
    rad = (int)(e[i]/1.414*12.0);
    if (mark[i] != 0) {
      Hv_DrawLine(xs-rad, ys-rad, xs+rad, ys+rad, Hv_black);
      Hv_DrawLine(xs+rad, ys-rad, xs-rad, ys+rad, Hv_black);
    }
  }    
  Hv_SetLineWidth(1);

  /* - draw a triangle at the impacts for each generated photon - */
  
  if (noMC) {
    for (i=0; i < noMC; i++){
      Hv_Point ximp;
      ximp.x = xc - (int)(xMC[i]*ppc);
      ximp.y = yc - (int)(yMC[i]*ppc);
      Hv_DrawSymbolsOnPoints(&ximp,1,10,Hv_blue,Hv_TRIANGLEUPSYMBOL);
    }
  }

  /*  checks run number and displays proper experimental setup */

  if (run < 5800) {
    Recoil(RecoilItem,region);
  }
  else {
    BSD(BSDItem,region);
    BGV(BGVItem,region);
  } 

  if (matchup == True) DrawMatchup(Item,region); /* draws particle match up if req. */
  if (lgscale == True) DrawLGscale(Item,region); /* draws scale for LGBs if req. */
  if (phodata == True) DrawConnect(Item,region); /* connects photons with their data */
}

/* - draw z axis of BSD on main view - */

void BSD(Hv_Item Item, Hv_Region region)
{
  int  icolor = Hv_white;
  int  i,j;

  for (i=0; i<24; i++) {
    int tdcon = 0;
    float angs=(float)i*0.2618+0.2618;
    float ange=(float)(i+1)*0.2618+0.2618;
    int xs=xc-(short)(190.0*cos(angs));
    int ys=yc-(short)(190.0*sin(angs));
    int xe=xc-(short)(190.0*cos(ange));
    int ye=yc-(short)(190.0*sin(ange));
    Hv_SetLineWidth(2);
    Hv_DrawLine(xs, ys, xe, ye, Hv_black);

    /* loop over all 24 sectors of the BSD */
    
    Hv_SetLineWidth(5);
    xs=xc-(short)(195.0*cos(angs));
    ys=yc-(short)(195.0*sin(angs));
    xe=xc-(short)(195.0*cos(ange));
    ye=yc-(short)(195.0*sin(ange));
    if (bsd_s[i] != 0){  

    /* draw straight scintillators of BSD, color-coded ADC values */

      float deltaE = bsd_s[i];
      if (calibrated) {
        if (deltaE > lb_bsdE) {
          icolor = deltaE*(float)(Hv_red - Hv_white)/ub_bsdE + Hv_white;
          if (icolor > Hv_white) icolor= Hv_white;
          if (icolor < Hv_red)  icolor= Hv_red;
          Hv_DrawLine(xs, ys, xe, ye, icolor);
        }
        for (j=0; j<bsd_s_times[i].multicount; j++) {
          float time = bsd_s_times[i].time[j];
	  if ((time > lb_bsdT) && (time < ub_bsdT)) tdcon = 1;
        }
      }
      else {
        if (bsd_s[i] > LB_ADC) {
          icolor = bsd_s[i]*(float)(Hv_red - Hv_white)/BSD_RED_HOT + Hv_white;
          if (icolor > Hv_white) icolor= Hv_white;
          if (icolor < Hv_red)  icolor= Hv_red;
          Hv_DrawLine(xs, ys, xe, ye, icolor);
        }
        for (j=0; j<bsd_s_times[i].multicount; j++) {
          int tdc = bsd_s_times[i].time[j];
	  if ((tdc > LB_BSDTDC) && (tdc < UB_BSDTDC)) tdcon = 1;
        }
      }
    }
    if (tdcon) {
      int xa=xs+(short)((xe-xs)/3.0);
      int ya=ys+(short)((ye-ys)/3.0);
      int xb=xs+(short)(2.0*(xe-xs)/3.0);
      int yb=ys+(short)(2.0*(ye-ys)/3.0);
      int xg=xs+(short)(2.0*(xe-xs)/5.0);
      int xh=xs+(short)(3.0*(xe-xs)/5.0);
      int yg=ys+(short)(2.0*(ye-ys)/5.0);
      int yh=ys+(short)(3.0*(ye-ys)/5.0);
      Hv_SetLineWidth(2);
      Hv_DrawLine(xa, ya, xb, yb, Hv_black);
      Hv_DrawLine(xg, yg, xh, yh, Hv_white);
    }
  }
  Hv_SetLineWidth(1);
}  /* end BSD */

/* draw the BGV upstream and downstream counters on main view */

void BGV(Hv_Item Item, Hv_Region region)
{
  int        icolor = Hv_white;
  int        i,j;

  for (i=0; i<24; i++) {
    int tdcon;
    float angs=(float)i*0.2618 - 0.1309;
    float ange=(float)(i+1)*0.2618 - 0.1309;
    int xs=xc-(short)(205.0*cos(angs));
    int ys=yc-(short)(205.0*sin(angs));
    int xe=xc-(short)(205.0*cos(ange));
    int ye=yc-(short)(205.0*sin(ange));

    Hv_SetLineWidth(2);
    Hv_DrawLine(xs, ys, xe, ye, Hv_black);

    /* loop over all 24 sectors of BGV */
    
    Hv_SetLineWidth(5);
    xs=xc-(short)(210.0*cos(angs));
    ys=yc-(short)(210.0*sin(angs));
    xe=xc-(short)(210.0*cos(ange));
    ye=yc-(short)(210.0*sin(ange));
    tdcon = 0;
    if (bgv_ds[i] != 0) {

 /* draw downstream BGV signals, color-coded ADC values */

      float energy = bgv_ds[i];
      if (calibrated) {
        if (energy > lb_bgvE) {
          icolor = energy*(float)(Hv_red - Hv_white)/ub_bgvE + Hv_white;
          if (icolor > Hv_white) icolor = Hv_white;
          if (icolor < Hv_red) icolor = Hv_red;
          Hv_DrawLine(xs, ys, xe, ye, icolor);
        }
        for (j=0; j<bgv_ds_times[i].multicount; j++) {
          float time = bgv_ds_times[i].time[j];
	  if ((time > lb_bgvT) && (time < ub_bgvT)) tdcon = 1;
        }
      }
      else {
        if (bgv_ds[i] > LB_ADC) {
          icolor = bgv_ds[i]*(float)(Hv_red - Hv_white)/BGV_RED_HOT + Hv_white;
          if (icolor > Hv_white) icolor = Hv_white;
          if (icolor < Hv_red) icolor = Hv_red;
          Hv_DrawLine(xs, ys, xe, ye, icolor);
        }
        for (j=0; j<bgv_ds_times[i].multicount; j++) {
          int tdc = bgv_ds_times[i].time[j];
	  if ((tdc > LB_BGVTDC) && (tdc < UB_BGVTDC)) tdcon = 1;
        }
      }
    }

    if (tdcon) {
      int xa=xs+(short)((xe-xs)/3.0);
      int ya=ys+(short)((ye-ys)/3.0);
      int xb=xs+(short)(2.0*(xe-xs)/3.0);
      int yb=ys+(short)(2.0*(ye-ys)/3.0);
      int xg=xs+(short)(2.0*(xe-xs)/5.0);
      int xh=xs+(short)(3.0*(xe-xs)/5.0);
      int yg=ys+(short)(2.0*(ye-ys)/5.0);
      int yh=ys+(short)(3.0*(ye-ys)/5.0);
      Hv_SetLineWidth(2);
      Hv_DrawLine(xa, ya, xb, yb, Hv_black);
      Hv_DrawLine(xg, yg, xh, yh, Hv_white);
    }

    Hv_SetLineWidth(5);
    xs=xc-(short)(220.0*cos(angs));
    ys=yc-(short)(220.0*sin(angs));
    xe=xc-(short)(220.0*cos(ange));
    ye=yc-(short)(220.0*sin(ange));
    tdcon = 0;
    if (bgv_us[i] != 0) {

 /* draw upstream BGV signals, color-coded ADC values */

      float energy = bgv_us[i];
      if (calibrated) {
        if (energy > lb_bgvE) {
          icolor = energy*(float)(Hv_red - Hv_white)/ub_bgvE + Hv_white;
          if (icolor > Hv_white) icolor = Hv_white;
          if (icolor < Hv_red) icolor = Hv_red;
          Hv_DrawLine(xs, ys, xe, ye, icolor);
        }
        for (j=0; j<bgv_us_times[i].multicount; j++) {
          float time = bgv_us_times[i].time[j];
	  if ((time > lb_bgvT) && (time < ub_bgvT)) tdcon = 1;
        }
      }
      else {
        if (bgv_us[i] > LB_ADC) {
          icolor = bgv_us[i]*(float)(Hv_red - Hv_white)/BGV_RED_HOT + Hv_white;
          if (icolor > Hv_white) icolor = Hv_white;
          if (icolor < Hv_red) icolor = Hv_red;
          Hv_DrawLine(xs, ys, xe, ye, icolor);
        }
        for (j=0; j<bgv_us_times[i].multicount; j++) {
          int tdc = bgv_us_times[i].time[j];
	  if ((tdc > LB_BGVTDC) && (tdc < UB_BGVTDC)) tdcon = 1;
        }
      }
    }

    if (tdcon) {
      int xa=xs+(short)((xe-xs)/3.0);
      int ya=ys+(short)((ye-ys)/3.0);
      int xb=xs+(short)(2.0*(xe-xs)/3.0);
      int yb=ys+(short)(2.0*(ye-ys)/3.0);
      int xg=xs+(short)(2.0*(xe-xs)/5.0);
      int xh=xs+(short)(3.0*(xe-xs)/5.0);
      int yg=ys+(short)(2.0*(ye-ys)/5.0);
      int yh=ys+(short)(3.0*(ye-ys)/5.0);
      Hv_SetLineWidth(2);
      Hv_DrawLine(xa, ya, xb, yb, Hv_black);
      Hv_DrawLine(xg, yg, xh, yh, Hv_white);
    }
  }
  Hv_SetLineWidth(1);
}  /* end BGV() */


/* - draw the proton recoil detector - */

void Recoil(Hv_Item Item, Hv_Region region)
{
  int		i, j;
  char		text[40];
  Hv_String	*Hvtext;
  
  /* loop over six sectors of hexagon  */

  for (i=0; i<6; i++) {
    float angs=(float)i*1.0472 + 0.2618;
    float ange=(float)(i+1)*1.0472 + 0.2618;
    int xs=xc-(short)(210.0*cos(angs));
    int ys=yc-(short)(210.0*sin(angs));
    int xe=xc-(short)(210.0*cos(ange));
    int ye=yc-(short)(210.0*sin(ange));
    Hv_SetLineWidth(2);
    Hv_DrawLine(xs, ys, xe, ye, Hv_black);

    /*  considerably hacked into here    dsa */   
    /* loop over all 12 sectors of RPD  */

    Hv_SetLineWidth(5);
    xs=xc-(short)(215.0*cos(angs));
    ys=yc-(short)(215.0*sin(angs));
    xe=xc-(short)(215.0*cos(ange));
    ye=yc-(short)(215.0*sin(ange));
    if (rpd_g[i] != 0)          /*  draw the g-counters, color-coded ADC values */
       Hv_DrawLine(xs, ys, xe, ye, Hv_white - rpd_g[i]/143);
    xs=xc-(short)(225.0*cos(angs));
    ys=yc-(short)(225.0*sin(angs));
    xe=xc-(short)(225.0*cos(ange));
    ye=yc-(short)(225.0*sin(ange));
    j = i*2;                    /* the even-numbered e-counter  */
    if (rpd_e[j] != 0)  
       Hv_DrawLine(xs, ys, (xs+xe)/2, (ys+ye)/2, Hv_white - rpd_e[j]/143);
    j = i*2 + 1;                /* the odd-numbered e-counter  */
    if (rpd_e[j] != 0)      
       Hv_DrawLine((xs+xe)/2, (ys+ye)/2, xe, ye, Hv_white - rpd_e[j]/143);

/*   this next deleted DSA
          sprintf(text, "Proton angle: %d-%d", (j  -1)*30, j  *30);
      Hvtext = Hv_CreateString(text);		   
      if (j   < 7)   
        Hv_StringDraw(xc-280+j  *60, yc-215, Hvtext, NULL);
      else
        Hv_StringDraw(xc-280+(13-j  )*60, yc+200, Hvtext, NULL);
      Hv_DestroyString(Hvtext);
*/
  }
  Hv_SetLineWidth(1);
}

void CPV(Hv_Item Item, Hv_Region region)
{
  char text[50];
  Hv_Rect area, rect, border;
  int bgcolor = Hv_powderBlue;
  int i,j;
  float CPV_w[3] = {17.188, 20.562, 22.188};
  float CPV_h[4]= {1.25, 3.5, 6.75,2.25};
  float	CPV_X[3]= {17.188,22.188,1.626};
  float	CPV_Y[8]= {22.625, 15.875, 12.375, 8.875, 5.375, 3.125,1.875,.625};

  /*in order of actual paddle nos */
  int   pos[30][4] = {{0,1,0,2},{0,1,0,2},{0,2,2,1},{1,2,2,1},{0,3,2,1},
		      {1,3,2,1},{0,4,2,1},{1,4,2,1},{0,5,2,3},{1,5,2,3},
		      {0,6,2,0},{1,6,2,0},{2,7,1,0},{1,7,1,0},{2,7,1,0},
		      {1,7,1,0},{2,6,1,0},{1,6,1,0},{1,5,2,0},{1,5,2,0},
		      {0,4,2,3},{1,4,2,3},{0,3,2,1},{1,3,2,1},{0,2,2,1},
		      {1,2,2,1},{0,1,2,1},{1,1,2,1},{0,0,0,2},{0,0,0,2}};

  float         ppi = 7.9558;  /*pixel per inch*/
  Hv_String	*Hvtext;

  Hv_CopyRect(&area, Hv_GetItemArea(Item));	/* gets the area */
  Hv_FillRect(&area, Hv_black);		/* clears the area */
  Hv_GetItemCenter(Item, &xctcpv, &yctcpv);

  Hv_SetRect(&border, xctcpv-200,yctcpv-215,400,35.000025);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border, xctcpv-200,yctcpv-179.999975,63.2557096,53.70165);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border,xctcpv+136.7442904,yctcpv-179.999975,63.2557096,53.70165);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border, xctcpv-200,yctcpv-126.298325,23.4767096,254.59665);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border,xctcpv+176.5232904,yctcpv-126.298325,23.4767096,254.59665);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border, xctcpv-200,yctcpv+179.999975,400,35.000025);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border, xctcpv-200,yctcpv+126.298325,63.2557096,53.70165);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border,xctcpv+136.7442904,yctcpv+126.298325,63.2557096,53.70165);
  Hv_FillRect(&border,bgcolor);
  Hv_SetRect(&border,xctcpv-12.110495,yctcpv-14.110495,29.22099,29.22099);
  Hv_FillRect(&border,bgcolor);
  sprintf(text, "Event number: %d", event);
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left+30, area.top+10, Hvtext, NULL);
 
  
  DrawCPV(xctcpv-(CPV_X[0]*ppi),yctcpv-(CPV_Y[0]*ppi),CPV_w[0]*ppi,CPV_h[2]*ppi,Hv_black);   /*29  top left "fat" paddle*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv-(CPV_Y[1]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);   /*27 down the left side*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv-(CPV_Y[2]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);   /*25*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv-(CPV_Y[3]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);   /*23*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv-(CPV_Y[4]*ppi),CPV_w[2]*ppi,CPV_h[3]*ppi,Hv_black);   /*21*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv-(CPV_Y[5]*ppi),CPV_w[2]*ppi,CPV_h[0]*ppi,Hv_black);   /*19*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv-(CPV_Y[6]*ppi),CPV_w[1]*ppi,CPV_h[0]*ppi,Hv_black);   /*17*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv-(CPV_Y[7]*ppi),CPV_w[1]*ppi,CPV_h[0]*ppi,Hv_black);   /*15*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv+(CPV_Y[7]*ppi),CPV_w[1]*ppi,CPV_h[0]*ppi,Hv_black);   /*13*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv+(CPV_Y[6]*ppi),CPV_w[2]*ppi,CPV_h[0]*ppi,Hv_black);   /*11*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv+(CPV_Y[5]*ppi),CPV_w[2]*ppi,CPV_h[3]*ppi,Hv_black);   /*9*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv+(CPV_Y[4]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);   /*7*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv+(CPV_Y[3]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);   /*5*/
  DrawCPV(xctcpv-(CPV_X[1]*ppi),yctcpv+(CPV_Y[2]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);   /*3*/
  DrawCPV(xctcpv-(CPV_X[0]*ppi),yctcpv+(CPV_Y[1]*ppi),CPV_w[0]*ppi,CPV_h[2]*ppi,Hv_black);   /*bottom left "fat" paddle 1 */
  DrawCPV(xctcpv,yctcpv+(CPV_Y[1]*ppi),CPV_w[0]*ppi,CPV_h[2]*ppi,Hv_black);                  /*bottom right "fat" paddle 0 */
  DrawCPV(xctcpv,yctcpv+(CPV_Y[2]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);                  /*up the right side 2 */
  DrawCPV(xctcpv,yctcpv+(CPV_Y[3]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);                  /*4*/
  DrawCPV(xctcpv,yctcpv+(CPV_Y[4]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);                  /*6*/
  DrawCPV(xctcpv,yctcpv+(CPV_Y[5]*ppi),CPV_w[2]*ppi,CPV_h[3]*ppi,Hv_black);                  /*8*/
  DrawCPV(xctcpv,yctcpv+(CPV_Y[6]*ppi),CPV_w[2]*ppi,CPV_h[0]*ppi,Hv_black);                  /*10*/
  DrawCPV(xctcpv+(CPV_X[2]*ppi),yctcpv+(CPV_Y[7]*ppi),CPV_w[1]*ppi,CPV_h[0]*ppi,Hv_black);   /*12*/
  DrawCPV(xctcpv+(CPV_X[2]*ppi),yctcpv-(CPV_Y[7]*ppi),CPV_w[1]*ppi,CPV_h[0]*ppi,Hv_black);   /*14*/
  DrawCPV(xctcpv+(CPV_X[2]*ppi),yctcpv-(CPV_Y[6]*ppi),CPV_w[1]*ppi,CPV_h[0]*ppi,Hv_black);   /*16*/
  DrawCPV(xctcpv,yctcpv-(CPV_Y[5]*ppi),CPV_w[2]*ppi,CPV_h[0]*ppi,Hv_black);                  /*18*/
  DrawCPV(xctcpv,yctcpv-(CPV_Y[4]*ppi),CPV_w[2]*ppi,CPV_h[3]*ppi,Hv_black);                  /*20*/
  DrawCPV(xctcpv,yctcpv-(CPV_Y[3]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);                  /*22*/
  DrawCPV(xctcpv,yctcpv-(CPV_Y[2]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);                  /*24*/
  DrawCPV(xctcpv,yctcpv-(CPV_Y[1]*ppi),CPV_w[2]*ppi,CPV_h[1]*ppi,Hv_black);                  /*26*/
  DrawCPV(xctcpv,yctcpv-(CPV_Y[0]*ppi),CPV_w[0]*ppi,CPV_h[2]*ppi,Hv_black);                  /*28 upper right paddle*/

  /* fill in CPV blocks with ADC value from color scale */

  for (i = 0; i < 30; i++){
    int icolor = Hv_white; 
    int tdcon = 0;
    if (cpv[i] > 0) {
      if (calibrated) {
        float deltaE = cpv[i];
        if (deltaE > lb_cpvE) {
          icolor = deltaE*(float)(Hv_red - Hv_white)/ub_cpvE + Hv_white;
	  if (icolor > Hv_white) icolor= Hv_white;
	  if (icolor < Hv_red) icolor= Hv_red;
        }
	for (j = 0; j < cpv_times[i].multicount; j++) {
          float time = cpv_times[i].time[j];
	  if ((time > lb_cpvT) && (time < ub_cpvT)) tdcon = 1;
        }
      }
      else {
        int adc = cpv[i];
        if (adc > LB_ADC) {
          icolor = adc*(float)(Hv_red - Hv_white)/CPV_RED_HOT + Hv_white;
	  if (icolor > Hv_white) icolor= Hv_white;
	  if (icolor < Hv_red) icolor= Hv_red;
        }
	for (j = 0; j < cpv_times[i].multicount; j++) {
          float tdc = cpv_times[i].time[j];
	  if ((tdc > LB_CPVTDC) && (tdc < UB_CPVTDC)) tdcon = 1;
        }
      }
    }

    switch (i) {  
    case 29:
      Hv_SetRect(&rect,(xctcpv-(CPV_X[pos[i][0]]*ppi)+1),
		 (yctcpv-(CPV_Y[pos[i][1]]*ppi)+1),
		 (CPV_w[pos[i][2]]*ppi)-1,(CPV_h[pos[i][3]]*ppi)-1);
      break;
    case 15:
      Hv_SetRect(&rect,(xctcpv-(CPV_X[pos[i][0]]*ppi)+1),
		 (yctcpv-(CPV_Y[pos[i][1]]*ppi)+1),
		 (CPV_w[pos[i][2]]*ppi),(CPV_h[pos[i][3]]*ppi)-1);
      break;
    case 17:
      Hv_SetRect(&rect,(xctcpv-(CPV_X[pos[i][0]]*ppi)+1),
		 (yctcpv-(CPV_Y[pos[i][1]]*ppi)+1),
		 (CPV_w[pos[i][2]]*ppi),(CPV_h[pos[i][3]]*ppi));
      break;
    case 19:
    case 21:
    case 23:
    case 25:
    case 27:
      Hv_SetRect(&rect,(xctcpv-(CPV_X[pos[i][0]]*ppi)+1),
		 (yctcpv-(CPV_Y[pos[i][1]]*ppi)+1),
		 (CPV_w[pos[i][2]]*ppi)-1,(CPV_h[pos[i][3]]*ppi));
      break;
    case 13:
      Hv_SetRect(&rect,(xctcpv-(CPV_X[pos[i][0]]*ppi)+1),
		 (yctcpv+(CPV_Y[pos[i][1]]*ppi)+1),
		 (CPV_w[pos[i][2]]*ppi),(CPV_h[pos[i][3]]*ppi));
      break;
    case 1: 
    case 3:
    case 5:
    case 7:
    case 9:
    case 11:
      Hv_SetRect(&rect,(xctcpv-(CPV_X[pos[i][0]]*ppi)+1),
		 (yctcpv+(CPV_Y[pos[i][1]]*ppi)+1),
		 (CPV_w[pos[i][2]]*ppi)-1,(CPV_h[pos[i][3]]*ppi));
      break;
    case 0:
    case 2:
    case 4:
    case 6:
    case 8:
    case 10:
      Hv_SetRect(&rect,xctcpv+1,(yctcpv+(CPV_Y[pos[i][1]]*ppi)+1),
		 (CPV_w[pos[i][2]]*ppi)-1,(CPV_h[pos[i][3]]*ppi));
      break;
    case 12:
      Hv_SetRect(&rect,xctcpv+(CPV_X[pos[i][0]]*ppi)+1,
		 yctcpv+(CPV_Y[pos[i][1]]*ppi)+1,(CPV_w[pos[i][2]]*ppi),
		 (CPV_h[pos[i][3]]*ppi));
      break;
    case 14:
      Hv_SetRect(&rect,(xctcpv+(CPV_X[pos[i][0]]*ppi))+1, 
		 (yctcpv-(CPV_Y[pos[i][1]]*ppi))+1,
		 (CPV_w[pos[i][2]]*ppi),(CPV_h[pos[i][3]]*ppi)-1);
      break;
    case 16:
      Hv_SetRect(&rect,(xctcpv+(CPV_X[pos[i][0]]*ppi))+1, 
		 (yctcpv-(CPV_Y[pos[i][1]]*ppi))+1,
		 (CPV_w[pos[i][2]]*ppi),(CPV_h[pos[i][3]]*ppi));
      break;
    case 28:
      Hv_SetRect(&rect,xctcpv+1,yctcpv-(CPV_Y[pos[i][1]]*ppi)+1,
		 (CPV_w[pos[i][2]]*ppi)-1,(CPV_h[pos[i][3]]*ppi)-1);
      break;
    case 18:
    case 20:
    case 22:
    case 24:
    case 26:
      Hv_SetRect(&rect,xctcpv+1,yctcpv-(CPV_Y[pos[i][1]]*ppi)+1,
		 (CPV_w[pos[i][2]]*ppi)-1,(CPV_h[pos[i][3]]*ppi));
      break;
    }

    if (tdcon == 1) {
      Hv_FillPatternRect(&rect, icolor, Hv_VSTRIPEPAT);
    }
    else {
      Hv_FillRect(&rect,icolor);
    }
  }
}

void DrawCPV(float Xstart,float Ystart,float width,float height,short icolor)
{
  Hv_DrawLine(Xstart,Ystart, Xstart+width,Ystart,icolor);
  Hv_DrawLine(Xstart+width,Ystart,Xstart+width,Ystart+height,icolor);
  Hv_DrawLine(Xstart+width,Ystart+height,Xstart,Ystart+height,icolor);
  Hv_DrawLine(Xstart,Ystart+height,Xstart,Ystart,icolor);
}

/* - draw the BSD in 2D (z,phi)- */

void BSDUnfolded(Hv_Item Item, Hv_Region region)
{
  Hv_Rect        area, rect;
  Hv_String      *Hvtext;

  char text[40];
  int  i,j;
 
  /*  horizontal distance to the border from the center */
  int     xbsd;
  /* horizontal and vertical length between the left/right counters  */
  int     LW, LH;
  /* horizontal distance between straight counters */
  int     SW;

  xbsd = 288;
  LW   = 48;
  LH   = 48;
  SW   = 24;
   
  Hv_CopyRect(&area, Hv_GetItemArea(Item));  /* gets the area */
  Hv_FillRect(&area, Hv_gray14);             /* clears the area */
  Hv_GetItemCenter(Item, &xctbsd, &yctbsd);
    
  sprintf(text, "Event number: %d", event);
  Hvtext = Hv_CreateString(text);
  Hv_StringDraw(area.left+30, area.top+10, Hvtext, NULL);
 
  /* fill in straight BSD counters with ADC value from color scale */

  for (i = 0; i < 24; i++) {
    int tdcon = 0;
    int icolor = Hv_white;
    if (calibrated) {
      float deltaE = bsd_s[i];
      if (deltaE > lb_bsdE) {
        icolor = deltaE*(float)(Hv_red - Hv_white)/ub_bsdE + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < bsd_s_times[i].multicount; j++) {
        float time = bsd_s_times[i].time[j];
        if ((time > lb_bsdT) && (time < ub_bsdT)) tdcon = 1;
      }
    }
    else {
      int adc = bsd_s[i];
      if (adc > LB_ADC) {
        icolor = adc*(float)(Hv_red - Hv_white)/BSD_RED_HOT + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < bsd_s_times[i].multicount; j++) {
        int tdc = bsd_s_times[i].time[j];
        if ((tdc > LB_BSDTDC) && (tdc < UB_BSDTDC)) tdcon = 1;
      }
    }
    Hv_SetRect(&rect,xctbsd-288+(SW*i),yctbsd-96,SW,192);
    if (tdcon) {
      Hv_FillPatternRect(&rect,icolor,Hv_HSTRIPEPAT);
    }
    else if (icolor) {
      Hv_FillRect(&rect,icolor);
    }
  }

  /* fill left counters with ADC value from color scale */

  for (i = 0; i < 12; i++){
    int tdcon = 0;
    int icolor = 0;
    if (calibrated) {
      float deltaE = bsd_l[i];
      if (deltaE > lb_bsdE) {
        icolor = deltaE*(float)(Hv_red - Hv_white)/ub_bsdE + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < bsd_l_times[i].multicount; j++) {
        float time = bsd_l_times[i].time[j];
        if ((time > lb_bsdT) && (time < ub_bsdT)) tdcon = 1;
      }
    }
    else {
      int adc = bsd_l[i];
      if (adc > LB_ADC) {
        icolor = adc*(float)(Hv_red - Hv_white)/BSD_RED_HOT + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < bsd_l_times[i].multicount; j++) {
        int tdc = bsd_l_times[i].time[j];
        if ((tdc > LB_BSDTDC) && (tdc < UB_BSDTDC)) tdcon = 1;
      }
    }

    switch(i) {
    case 0:
    case 1:
    case 2:
    case 3:
      if (tdcon)  {
        PatternFill(xctbsd-288-LW,yctbsd+LH-(LH*i),1,Hv_black,icolor);
        PatternFill(xctbsd+96+(LW*i),yctbsd-96,1,Hv_black,icolor);
      }
      else if (icolor) {
        PatternFill(xctbsd-288-LW,yctbsd+LH-(LH*i),1,icolor,icolor);
        PatternFill(xctbsd+96+(LW*i),yctbsd-96,1,icolor,icolor);
      }
      break;
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
      if (tdcon) {
        PatternFill(xctbsd-288+LW*(i-4),yctbsd-96,1,Hv_black,icolor);
      }
      else if (icolor) {
        PatternFill(xctbsd-288+LW*(i-4),yctbsd-96,1,icolor,icolor);
      }
      break;
    }
  }

  /* fill right counters with ADC value from color scale */

  for (i = 0; i < 12; i++) {
    int tdcon = 0;
    int icolor = 0;
    if (calibrated) {
      float deltaE = bsd_r[i];
      if (deltaE > lb_bsdE) {
        icolor = deltaE*(float)(Hv_red - Hv_white)/ub_bsdE + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < bsd_r_times[i].multicount; j++) {
        float time = bsd_r_times[i].time[j];
        if ((time > lb_bsdT) && (time < ub_bsdT)) tdcon = 1;
      }
    }
    else {
      int adc = bsd_r[i];
      if (adc > LB_ADC) {
        icolor = adc*(float)(Hv_red - Hv_white)/BSD_RED_HOT + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < bsd_r_times[i].multicount; j++) {
        int tdc = bsd_r_times[i].time[j];
        if ((tdc > LB_BSDTDC) && (tdc < UB_BSDTDC)) tdcon = 1;
      }
    }

    switch (i) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      if (tdcon) {
        PatternFill(xctbsd-96+(LW*i),yctbsd-96,-1,Hv_black,icolor);
      }
      else if (icolor) {
        PatternFill(xctbsd-96+(LW*i),yctbsd-96,-1,icolor,icolor);
      }
      break;
    case 8:
    case 9:
    case 10:
    case 11:
      if (tdcon) {
        PatternFill(xctbsd+288,yctbsd+LH*(i-10),-1,Hv_black,icolor);
        PatternFill(xctbsd-288+LW*(i-8), yctbsd-96,-1,Hv_black,icolor);
      }
      else if (icolor) {
        PatternFill(xctbsd+288,yctbsd+LH*(i-10),-1,icolor,icolor);
        PatternFill(xctbsd-288+LW*(i-8), yctbsd-96,-1,icolor,icolor);
      }
      break;
    }
  }

  /* fill in straight BSD counters with hatches for a valid TDC value */

  for (i = 0; i < 24; i++) {
    int tdcon = 0;
    int icolor = Hv_white;
    if (calibrated) {
      float deltaE = bsd_s[i];
      if (deltaE > lb_bsdE) {
        icolor = deltaE*(float)(Hv_red - Hv_white)/ub_bsdE + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
        for (j = 0; j < bsd_s_times[i].multicount; j++) {
          float time = bsd_s_times[i].time[j];
          if ((time > lb_bsdT) && (time < ub_bsdT)) tdcon = 1;
        }
      }
    }
    else {
      int adc = bsd_s[i];
      if (adc > LB_ADC) {
        icolor = adc*(float)(Hv_red - Hv_white)/BSD_RED_HOT + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
        for (j = 0; j < bsd_s_times[i].multicount; j++) {
          int tdc = bsd_s_times[i].time[j];
          if ((tdc > LB_BSDTDC) && (tdc < UB_BSDTDC)) tdcon = 1;
        }
      }
    }
    Hv_SetRect(&rect,xctbsd-288+(SW*i),yctbsd-96,SW,192);
    if (tdcon) {
      Hv_FillPatternRect(&rect,icolor,Hv_HSTRIPEPAT);
    }
  }

/* fill left counters with hatching for TDC*/

  for (i = 0; i < 12; i++){
    int tdcon = 0;
    int icolor = 0;
    if (calibrated) {
      float deltaE = bsd_l[i];
      if (deltaE > lb_bsdE) {
        icolor = deltaE*(float)(Hv_red - Hv_white)/ub_bsdE + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
        for (j = 0; j < bsd_l_times[i].multicount; j++) {
          float time = bsd_l_times[i].time[j];
          if ((time > lb_bsdT) && (time < ub_bsdT)) tdcon = 1;
        }
      }
    }
    else {
      int adc = bsd_l[i];
      if (adc > LB_ADC) {
        icolor = adc*(float)(Hv_red - Hv_white)/BSD_RED_HOT + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
        for (j = 0; j < bsd_l_times[i].multicount; j++) {
          int tdc = bsd_l_times[i].time[j];
          if ((tdc > LB_BSDTDC) && (tdc < UB_BSDTDC)) tdcon = 1;
        }
      }
    }

    switch(i) {
    case 0:
    case 1:
    case 2:
    case 3:
      if (tdcon)  {
        PatternFill(xctbsd-288-LW,yctbsd+LH-(LH*i),1,Hv_black,icolor);
        PatternFill(xctbsd+96+(LW*i),yctbsd-96,1,Hv_black,icolor);
      }
      break;
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
      if (tdcon) {
        PatternFill(xctbsd-288+LW*(i-4),yctbsd-96,1,Hv_black,icolor);
      }
      break;
    }
  }

/* fill right counters with hatching for TDC */

  for (i = 0; i < 12; i++) {
    int tdcon = 0;
    int icolor = 0;
    if (calibrated) {
      float deltaE = bsd_r[i];
      if (deltaE > lb_bsdE) {
        icolor = deltaE*(float)(Hv_red - Hv_white)/ub_bsdE + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
        for (j=0; j < bsd_r_times[i].multicount; j++) {
          float time = bsd_r_times[i].time[j];
          if ((time > lb_bsdT) && (time < ub_bsdT)) tdcon = 1;
        }
      }
    }
    else {
      int adc = bsd_r[i];
      if (adc > LB_ADC) {
        icolor = adc*(float)(Hv_red - Hv_white)/BSD_RED_HOT + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
        for (j = 0; j < bsd_r_times[i].multicount; j++) {
          int tdc = bsd_r_times[i].time[j];
          if ((tdc > LB_BSDTDC) && (tdc < UB_BSDTDC)) tdcon = 1;
        }
      }
    }

    switch (i) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      if (tdcon) {
        PatternFill(xctbsd-96+(LW*i),yctbsd-96,-1,Hv_black,icolor);
      }
      break;
    case 8:
    case 9:
    case 10:
    case 11:
      if (tdcon) {
        PatternFill(xctbsd+288,yctbsd+LH*(i-10),-1,Hv_black,icolor);
        PatternFill(xctbsd-288+LW*(i-8), yctbsd-96,-1,Hv_black,icolor);
      }
      break;
    }
  }

 /* Draw border around the BSD display */

  Hv_DrawLine(xctbsd-xbsd,yctbsd-96,xctbsd+xbsd,yctbsd-96,Hv_black);
  Hv_DrawLine(xctbsd+288,yctbsd-96,xctbsd+xbsd,yctbsd+96,Hv_black);
  Hv_DrawLine(xctbsd+288,yctbsd+96,xctbsd-xbsd,yctbsd+96,Hv_black);
  Hv_DrawLine(xctbsd-xbsd,yctbsd+96,xctbsd-xbsd,yctbsd-96,Hv_black);
  

 /* Draw the straight counters of the BSD */

  for (i = 0; i < 24; i++) {
    xbsd -= 24;  /* distance between straight counters */
    Hv_DrawLine(xctbsd-xbsd, yctbsd-96, xctbsd-xbsd, yctbsd+96, Hv_black);
  }

 /* Draw the left counters of the BSD */
 /* 48 is the vertical distance between the left counters */
 /* 48 is the horizontal distance between the left counters */
 
  for (i=0; i<15; i++) {
    int xs, ys, xe, ye;
    if (i <= 2) {
      xs = xctbsd-288;
      ys = yctbsd+LH-(LH*i);
      xe = xctbsd-288+LW*(i+1);
      ye = yctbsd+96;
    }
    else if (i <= 11) {
      xs = xctbsd-288+LW*(i-3);
      ys = yctbsd-96;
      xe = xctbsd-288+LW*(i+1);
      ye = yctbsd+96;
    }
    else {
      xs = xctbsd-288+LW*(i-3);
      ys = yctbsd-96;
      xe = xctbsd+288;
      ye = yctbsd+96-(LH*(i-11));
    }
    Hv_DrawLine(xs, ys, xe, ye, Hv_black);
  }
   
 /* Draw the right counters of the BSD */
 /* 48 is the vertical distance between the right counters */
 /* 48 is the horizontal distance between the right counters */
 
  for (i = 0; i < 15; i++) {
    int xs, ys, xe, ye;
    if (i <= 2) {
      xs = xctbsd-288;
      ys = yctbsd-LH+(LH*i);
      xe = xctbsd-288+LW*(i+1);
      ye = yctbsd-96;
    }
    else if ((i > 2) && (i <= 11)) {
      xs = xctbsd-288+LW*(i-3);
      ys = yctbsd+96;
      xe = xctbsd-288+LH*(i+1);
      ye = yctbsd-96;
    }
    else {
      xs = xctbsd+(LW*(i-9));
      ys = yctbsd+96;
      xe = xctbsd+288;
      ye = yctbsd-96+(LH*(i-11));
    }
    Hv_DrawLine(xs, ys, xe, ye, Hv_black);
  }

} /* end BSDUnfolded */


void PatternFill(int lx,int ty,int incx,short color_stripe,short color_reg)
{
/* function that will fill polygons with stripes for TDC in L or R counters */
/* and ADC value from color scale                                           */
/* Notes:                                                                   */
/* RIGHT = right most pixel                                                 */
/* BOTTOM = bottom most pixel                                               */
/* LW = the width of the parallelogram from left to right                   */
/* n = a shade is drawn every n lines                                       */

/* Parameters:                                                              */
/* ty = top most y in parallelogram                                         */
/* lx = what would be the left top most x of the parallelogram (may be < 0) */
/* color_stripe = color of the hatching                                     */
/* color_reg = color of the ADC color scale                                 */ 

  /* Note that if this is a left counter incx = 1
                           right counter incx = -1 */

  int x, y, sx, ex, corr_color;
  int RIGHT = xctbsd + 288;
  int BOTTOM = yctbsd + 96;
  int LW = 48;
  int n = 5;

  x = lx;
  for(y=ty; y < BOTTOM; y++) {
    if (((incx == 1) && !(y%n))||((incx == -1) && !((y+n/2)%n)))
      corr_color = color_stripe;
    else
      corr_color=color_reg;

    if(x < 0)
      sx = 0;
    else
      sx = x;

    if(x+LW > RIGHT)
      ex = RIGHT;
    else
      ex = x + LW;

    Hv_DrawLine(sx,y,ex,y,corr_color);
    x+=incx;
  } /* end for */
} /* end PatternFill */

void UPV(Hv_Item Item, Hv_Region region)
{
  Hv_Rect	area, rect, rect1, rect2, border;
   int i,j;
   int bgcolor = Hv_aliceBlue;
   char	 text[40];
   float UPV_w[2]= {7.1, 14.5};
   float UPV_h[3]= {6.74, 7.1, 14.5};
   float UPV_X[3]= {2.25, 3.55, 16.75};
   float UPV_Y[3]= {3.37, 10.47,16.75 };
   int   pos[6][4] = {{0,0,1,1},
		      {2,0,1,1},
		      {0,0,1,0},
		      {2,0,1,0},
		      {0,1,1,1},
		      {2,1,1,1}};
  
  float         ppi = 5.54;  /*pixel per inch*/
  Hv_String	*Hvtext;

  Hv_FillRect(&area, Hv_black);		/* clears the area */ 
  Hv_CopyRect(&area, Hv_GetItemArea(Item));	/* gets the area */
  Hv_GetItemCenter(Item, &xctupv, &yctupv);	

  Hv_SetRect(&border, xctupv-100, yctupv-125,200,32.205);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border, xctupv-100, yctupv-92.792, 80.333, 34.794);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border, xctupv+19.667, yctupv-92.792, 80.333, 34.794);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border,xctupv-100,yctupv-58.998,7.205,118.002);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border,xctupv+92.795,yctupv-58.998,7.205,118.002);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border,xctupv-100,yctupv+58.998,80.333,37.794);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border,xctupv+19.667,yctupv+58.998,80.333,37.794);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border, xctupv-100, yctupv+94.795,200,32.205);
  Hv_FillRect(&border, bgcolor);
  Hv_SetRect(&border, xctupv-12.465, yctupv-12.465,25.93,25.93);
  Hv_FillRect(&border, bgcolor);
  
  sprintf(text, "Event number: %d", event);
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left+30, area.top+10, Hvtext, NULL);
 
  
  /* the horizontal paddles are easy */
  DrawUPV(xctupv-(UPV_X[2]*ppi),yctupv-(UPV_Y[1]*ppi),UPV_w[1]*ppi,UPV_h[1]*ppi,Hv_black);   /*H5 horizontal down the left side*/
  DrawUPV(xctupv-(UPV_X[2]*ppi),yctupv-(UPV_Y[0]*ppi),UPV_w[1]*ppi,UPV_h[0]*ppi,Hv_black);   /*H3*/
  DrawUPV(xctupv-(UPV_X[2]*ppi),yctupv+(UPV_Y[0]*ppi),UPV_w[1]*ppi,UPV_h[1]*ppi,Hv_black);   /*H1*/
  DrawUPV(xctupv+(UPV_X[0]*ppi),yctupv+(UPV_Y[0]*ppi),UPV_w[1]*ppi,UPV_h[1]*ppi,Hv_black);   /*H0 horizontal up the right side*/
  DrawUPV(xctupv+(UPV_X[0]*ppi),yctupv-(UPV_Y[0]*ppi),UPV_w[1]*ppi,UPV_h[0]*ppi,Hv_black);   /*H2*/
  DrawUPV(xctupv+(UPV_X[0]*ppi),yctupv-(UPV_Y[1]*ppi),UPV_w[1]*ppi,UPV_h[1]*ppi,Hv_black);   /*H4*/
  
  /*time for the vertical paddles*/
  Hv_DrawLine((xctupv-(UPV_X[1])*ppi),(yctupv+(UPV_Y[2])*ppi),(xctupv+(UPV_X[1])*ppi),(yctupv+(UPV_Y[2])*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[1]*ppi)),(yctupv+(UPV_Y[2])*ppi),(xctupv+(UPV_X[1])*ppi),(yctupv+(10.65)*ppi),Hv_black);
  Hv_DrawLine((xctupv-(UPV_X[1]*ppi)),(yctupv+(UPV_Y[2])*ppi),(xctupv-(UPV_X[1])*ppi),(yctupv+(10.65)*ppi),Hv_black);
  Hv_DrawLine((xctupv-(UPV_X[1]*ppi)),(yctupv+(10.65)*ppi),(xctupv-(UPV_X[0])*ppi),(yctupv+(10.65)*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[1]*ppi)),(yctupv+(10.65)*ppi),(xctupv+(UPV_X[0])*ppi),(yctupv+(10.65)*ppi),Hv_black);
  Hv_DrawLine((xctupv-(UPV_X[0]*ppi)),(yctupv+(10.65)*ppi),(xctupv-(UPV_X[0])*ppi),(yctupv+(2.25)*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[0]*ppi)),(yctupv+(10.65)*ppi),(xctupv+(UPV_X[0])*ppi),(yctupv+(2.25)*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[0]*ppi)),(yctupv+(2.25)*ppi),(xctupv-(UPV_X[0]*ppi)),(yctupv+(2.25)*ppi),Hv_black);
  
  Hv_DrawLine((xctupv-(UPV_X[1])*ppi),(yctupv-(UPV_Y[2])*ppi),(xctupv+(UPV_X[1])*ppi),(yctupv-(UPV_Y[2])*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[1]*ppi)),(yctupv-(UPV_Y[2])*ppi),(xctupv+(UPV_X[1])*ppi),(yctupv-(9.65)*ppi),Hv_black);
  Hv_DrawLine((xctupv-(UPV_X[1]*ppi)),(yctupv-(UPV_Y[2])*ppi),(xctupv-(UPV_X[1])*ppi),(yctupv-(9.65)*ppi),Hv_black);
  Hv_DrawLine((xctupv-(UPV_X[1]*ppi)),(yctupv-(9.65)*ppi),(xctupv-(UPV_X[0])*ppi),(yctupv-(9.65)*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[1]*ppi)),(yctupv-(9.65*ppi)),(xctupv+(UPV_X[0]*ppi)),(yctupv-(9.65*ppi)),Hv_black);
  Hv_DrawLine((xctupv-(UPV_X[0]*ppi)),(yctupv-(9.65)*ppi),(xctupv-(UPV_X[0])*ppi),(yctupv-(2.25)*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[0]*ppi)),(yctupv-(9.65)*ppi),(xctupv+(UPV_X[0])*ppi),(yctupv-(2.25)*ppi),Hv_black);
  Hv_DrawLine((xctupv+(UPV_X[0]*ppi)),(yctupv-(2.25)*ppi),(xctupv-(UPV_X[0]*ppi)),(yctupv-(2.25)*ppi),Hv_black);

  
  /* fill in UPV blocks with ADC value from color scale */

  for (i = 0; i < 8; i++) {
    int tdcon = 0;
    int icolor = Hv_white;
    if (calibrated) {
      float deltaE = upv[i];
      if (deltaE > lb_upvE) {
        icolor = deltaE*(float)(Hv_red - Hv_white)/ub_upvE + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < upv_times[i].multicount; j++){
        float time = upv_times[i].time[j];
	if ((time > lb_upvT) && (time < ub_upvT)) tdcon = 1;
      }
    }
    else {
      int adc = upv[i];
      if (adc > LB_ADC) {
        icolor = adc*(float)(Hv_red - Hv_white)/UPV_RED_HOT + Hv_white;
        if (icolor > Hv_white) icolor= Hv_white;
        if (icolor < Hv_red) icolor= Hv_red;
      }
      for (j = 0; j < upv_times[i].multicount; j++) {
        int tdc = upv_times[i].time[j];
	if ((tdc > LB_UPVTDC) && (tdc < UB_UPVTDC)) tdcon = 1;
      }
    }

    switch (i) {  
    case 5:
      Hv_SetRect(&rect,(xctupv-(UPV_X[pos[i][0]]*ppi)+1),
		 (yctupv-(UPV_Y[pos[i][1]]*ppi)+1),
		 (UPV_w[pos[i][2]]*ppi)-1,(UPV_h[pos[i][3]]*ppi));
      if (tdcon) {
        Hv_FillPatternRect(&rect, icolor, Hv_VSTRIPEPAT);
      }
      else if (icolor) {
        Hv_FillRect(&rect,icolor);
      }    
      break;
    case 3:
      Hv_SetRect(&rect,(xctupv-(UPV_X[pos[i][0]]*ppi)+1),
		 (yctupv-(UPV_Y[pos[i][1]]*ppi)+1),
		 (UPV_w[pos[i][2]]*ppi)-1,(UPV_h[pos[i][3]]*ppi)-1);
      if (tdcon) {
        Hv_FillPatternRect(&rect, icolor, Hv_VSTRIPEPAT);
      }
      else if (icolor) {
        Hv_FillRect(&rect,icolor);
      }    
      break;
    case 2:
      Hv_SetRect(&rect,(xctupv+(UPV_X[pos[i][0]]*ppi)+1),
		 (yctupv-(UPV_Y[pos[i][1]]*ppi)+1),
		 (UPV_w[pos[i][2]]*ppi)-1,(UPV_h[pos[i][3]]*ppi)-1);
      if (tdcon) {
        Hv_FillPatternRect(&rect, icolor, Hv_VSTRIPEPAT);
      }
      else if (icolor) {
        Hv_FillRect(&rect,icolor);
      }    
      break;
    case 4:
      Hv_SetRect(&rect,(xctupv+(UPV_X[pos[i][0]]*ppi)+1),
		 (yctupv-(UPV_Y[pos[i][1]]*ppi)+1),
		 (UPV_w[pos[i][2]]*ppi)-1,(UPV_h[pos[i][3]]*ppi));
      if (tdcon) {
        Hv_FillPatternRect(&rect, icolor, Hv_VSTRIPEPAT);
      }
      else if (icolor) {
        Hv_FillRect(&rect,icolor);
      }    
      break;
    case 1:
      Hv_SetRect(&rect,(xctupv-(UPV_X[pos[i][0]]*ppi)+1),
		 (yctupv+(UPV_Y[pos[i][1]]*ppi)+1),
		 (UPV_w[pos[i][2]]*ppi)-1,(UPV_h[pos[i][3]]*ppi));
      if (tdcon) {
        Hv_FillPatternRect(&rect, icolor, Hv_VSTRIPEPAT);
      }
      else if (icolor) {
        Hv_FillRect(&rect,icolor);
      }    
      break;
    case 0:
       Hv_SetRect(&rect,(xctupv+(UPV_X[pos[i][0]]*ppi)+1),
		 (yctupv+(UPV_Y[pos[i][1]]*ppi)+1),
		 (UPV_w[pos[i][2]]*ppi)-1,(UPV_h[pos[i][3]]*ppi));
      if (tdcon) {
        Hv_FillPatternRect(&rect, icolor, Hv_VSTRIPEPAT);
      }
      else if (icolor) {
        Hv_FillRect(&rect,icolor);
      }    
      break;
    case 6:
      Hv_SetRect(&rect2,(xctupv-UPV_X[0]*ppi)+1,(yctupv+(2.25)*ppi)+1,(4.5*ppi),(8.22*ppi)+1);
      Hv_SetRect(&rect1,(xctupv-UPV_X[1]*ppi)+1,(yctupv+(10.65)*ppi),(UPV_w[0]*ppi)-1,(6.1*ppi)+1);
      if (tdcon) { 
	Hv_FillPatternRect(&rect1, icolor, Hv_VSTRIPEPAT);
	Hv_FillPatternRect(&rect2, icolor, Hv_VSTRIPEPAT);
      }
      else {
	Hv_FillRect(&rect1,icolor);
	Hv_FillRect(&rect2,icolor);
      }
      break;
    case 7:
      Hv_SetRect(&rect1,(xctupv-UPV_X[1]*ppi)+1,(yctupv-UPV_Y[2]*ppi)+1,(UPV_w[0]*ppi)-1,(6.1*ppi)-1);
      Hv_SetRect(&rect2,(xctupv-UPV_X[0]*ppi)+1,(yctupv-(10.65)*ppi),(4.5*ppi),(8.22*ppi)+2);
      if (tdcon){ 
	Hv_FillPatternRect(&rect1, icolor, Hv_VSTRIPEPAT);
	Hv_FillPatternRect(&rect2, icolor, Hv_VSTRIPEPAT);
      }
      else {
	Hv_FillRect(&rect1,icolor);
	Hv_FillRect(&rect2,icolor);
      }
      break;
    }
  }
}

void DrawUPV(float Xstart,float Ystart,float width,float height,short icolor)
{
  Hv_DrawLine(Xstart,Ystart, Xstart+width,Ystart,icolor);
  Hv_DrawLine(Xstart+width,Ystart,Xstart+width,Ystart+height,icolor);
  Hv_DrawLine(Xstart+width,Ystart+height,Xstart,Ystart+height,icolor);
  Hv_DrawLine(Xstart,Ystart+height,Xstart,Ystart,icolor);
}


/* - draw the photon tagger  - */

void Tagger(Hv_Item Item, Hv_Region region)
{
  Hv_Rect	area, rect;
  Hv_String	*Hvtext;

  int  bgcolor = Hv_white-5;
  char text[40];
  int  i;


  Hv_CopyRect(&area, Hv_GetItemArea(Item));	/* gets the area */
  Hv_FillRect(&area, bgcolor);			/* clears the area */
  Hv_GetItemCenter(Item, &xctag, &yctag);
 
  sprintf(text, "Event number: %d", event);
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left+30, area.top+10, Hvtext, NULL);

  xctag=xctag-nsize*7;

  /*   draw the tagger hits   */

  for (i=0; i<19; i++) {
    int k;
    int icolor = Hv_white;
    Hv_SetRect(&rect, xctag+nsize*i, yctag, nsize, nsize*2);
    for (k = 0; k < tagger_times[i].multicount; k++) {
      if (calibrated) {
        float ttag = tagger_times[i].time[k];
        if ((ttag > lb_tagT) && (ttag < ub_tagT)) {
          icolor = fabs(ttag)*(Hv_red - Hv_white)/100. + Hv_white;
          if (icolor > Hv_white) icolor = Hv_white;
          if (icolor < Hv_red) icolor = Hv_red;
        }
      }
      else {
        int ttag = tagger_times[i].time[k];
        if ((ttag < 400) && (ttag > 200)) {
          icolor = ttag*(float)(Hv_red - Hv_white)/TAGGER_FULL_SCALE + Hv_white;
          if (icolor > Hv_white) icolor = Hv_white;
          if (icolor < Hv_red) icolor = Hv_red;
        }
      }
    }
    Hv_FillRect(&rect, icolor);   /* scale: white - red */
  }
  for (i=19; i<38; i++) {
    int k;
    int icolor = Hv_white;
    Hv_SetRect(&rect, xctag+nsize*(i-19), yctag-nsize*2, nsize, nsize*2);
    for (k = 0; k < tagger_times[i].multicount; k++) {
      if (calibrated) {
        float ttag = tagger_times[i].time[k];
        if ((ttag > lb_tagT) && (ttag < ub_tagT)) {
          icolor = fabs(ttag)*(Hv_red - Hv_white)/100. + Hv_white;
          if (icolor > Hv_white) icolor= Hv_white;
          if (icolor < Hv_red) icolor= Hv_red;
        }
      }
      else {
        int ttag = tagger_times[i].time[k];
        if ((ttag < 320) && (ttag > 120)) {
          icolor = ttag*(float)(Hv_red - Hv_white)/TAGGER_FULL_SCALE + Hv_white;
          if (icolor > Hv_white) icolor = Hv_white;
          if (icolor < Hv_red) icolor = Hv_red;
        }
      }
    }
    Hv_FillRect(&rect, icolor);   /* scale: white - red */
  }

  /* draw the lines separating the tagger elements  */

  for (i=0; i<20; i++) {
     Hv_DrawLine(xctag+nsize*i, yctag-nsize*2, xctag+nsize*i, yctag+nsize*2, Hv_black);
     Hv_DrawLine(xctag+nsize*i, yctag-nsize*2, xctag+nsize*i, yctag+nsize*2, Hv_black);
  }

  Hv_DrawLine(xctag, yctag-nsize*2, xctag+nsize*19 , yctag-nsize*2, Hv_black);
  Hv_DrawLine(xctag, yctag, xctag+nsize*19, yctag, Hv_chocolate);
  Hv_DrawLine(xctag, yctag+nsize*2, xctag+nsize*19 , yctag+nsize*2, Hv_black);

  sprintf(text, "Left T-counters");
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(xctag+2*nsize, yctag+nsize*3, Hvtext, NULL);

  sprintf(text, "Right T-counters");
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(xctag+2*nsize, yctag-nsize*3.5, Hvtext, NULL);

  sprintf(text, "Low Energy");
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(xctag+nsize*19.5, yctag-nsize, Hvtext, NULL);

  sprintf(text, "High Energy");
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(xctag-nsize*6.5, yctag-nsize, Hvtext, NULL);
}

/* - draw particle match up - */

void DrawMatchup(Hv_Item Item, Hv_Region region)
{
  Hv_Rect	rect;
  int		i, al, at;
  char		text[40];
  Hv_String	*Hvtext;
  
  Hv_SetLineWidth(2);
  for (i=0; i < npi; i++) {
    Hv_SetRect(&rect, xc-190-i*20, yc-175+i*46, 105, 43);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);

    sprintf(text, "E: %4.3f GeV", epi[i]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+2, Hvtext, NULL);
    sprintf(text, "M: %4.3f GeV", mpi[i]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+15, Hvtext, NULL);
    sprintf(text, "t: %4.3f GeV\\s\\-\\^2", tpi[i]);
    Hvtext = Hv_CreateString(text);		   
    Hv_DrawCompoundString(rect.left+5, rect.top+28, Hvtext);
    Hv_DestroyString(Hvtext);
  
    Hv_DrawLine(rect.right, rect.bottom, xc-(int)(x[p1pi[i]]*ppc), yc-(int)(y[p1pi[i]]*ppc), Hv_black);
    Hv_DrawLine(rect.right, rect.bottom, xc-(int)(x[p2pi[i]]*ppc), yc-(int)(y[p2pi[i]]*ppc), Hv_black);
  
    Hv_SetRect(&rect, xc-205-i*20, yc-190+i*46, 30, 17);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);
    Hvtext = Hv_CreateString("Pi\\s\\-\\^0");		   
    Hv_DrawCompoundString(rect.left+4, rect.top+3, Hvtext);
    Hv_DestroyString(Hvtext);
  }
  if (neta != 0) {
    Hv_SetRect(&rect, xc-230, yc+86, 105, 43);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);

    sprintf(text, "E: %4.3f GeV", eta[0]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+2, Hvtext, NULL);
    sprintf(text, "M: %4.3f GeV", eta[1]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+15, Hvtext, NULL);
    sprintf(text, "t: %4.3f GeV\\s\\-\\^2", eta[2]);
    Hvtext = Hv_CreateString(text);		   
    Hv_DrawCompoundString(rect.left+5, rect.top+28, Hvtext);
    Hv_DestroyString(Hvtext);
    
    Hv_DrawLine(rect.right, rect.top, xc-(int)(x[p1eta]*ppc), yc-(int)(y[p1eta]*ppc), Hv_black);
    Hv_DrawLine(rect.right, rect.top, xc-(int)(x[p2eta]*ppc), yc-(int)(y[p2eta]*ppc), Hv_black);
  
    Hv_SetRect(&rect, xc-237, yc+127, 30, 19);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);
    Hvtext = Hv_CreateString("Eta");		   
    Hv_DrawCompoundString(rect.left+4, rect.top+3, Hvtext);
    Hv_DestroyString(Hvtext);
  }
  if (nomega != 0) {
    Hv_SetRect(&rect, xc-210, yc+132, 105, 43);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);

    sprintf(text, "E: %4.3f GeV", Omega[0]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+2, Hvtext, NULL);
    sprintf(text, "M: %4.3f GeV", Omega[1]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+15, Hvtext, NULL);
    sprintf(text, "t: %4.3f GeV\\s\\-\\^2", Omega[2]);
    Hvtext = Hv_CreateString(text);		   
    Hv_DrawCompoundString(rect.left+5, rect.top+28, Hvtext);
    Hv_DestroyString(Hvtext);
  
    Hv_DrawLine(rect.right, rect.top, xc-85, yc-132, Hv_black);
    Hv_DrawLine(rect.right, rect.top, xc-(int)(x[p2omega]*ppc), yc-(int)(y[p2omega]*ppc), Hv_black);
  
    Hv_SetRect(&rect, xc-217, yc+173, 42, 19);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);
    Hvtext = Hv_CreateString("Omega");		   
    Hv_DrawCompoundString(rect.left+3, rect.top+3, Hvtext);
    Hv_DestroyString(Hvtext);
  }
  if (nphi != 0) {
    Hv_SetRect(&rect, xc-210, yc+132, 105, 43);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);

    sprintf(text, "E: %4.3f GeV", phi[0]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+2, Hvtext, NULL);
    sprintf(text, "M: %4.3f GeV", phi[1]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(rect.left+5, rect.top+15, Hvtext, NULL);
    sprintf(text, "t: %4.3f GeV\\s\\-\\^2", phi[2]);
    Hvtext = Hv_CreateString(text);		   
    Hv_DrawCompoundString(rect.left+5, rect.top+28, Hvtext);
    Hv_DestroyString(Hvtext);
  
    Hv_DrawLine(rect.right, rect.top, xc-85, yc-132, Hv_black);
    Hv_DrawLine(rect.right, rect.top, xc-(int)(x[p2phi]*ppc), yc-(int)(y[p2phi]*ppc), Hv_black);
  
    Hv_SetRect(&rect, xc-217, yc+173, 30, 19);
    Hv_FillRect(&rect, Hv_white);
    Hv_FrameRect(&rect, Hv_black);
    Hvtext = Hv_CreateString("Phi");		   
    Hv_DrawCompoundString(rect.left+3, rect.top+3, Hvtext);
    Hv_DestroyString(Hvtext);
  }
  Hv_SetLineWidth(1);
}

/* - draw a scale for the LGBs - */

void DrawLGscale(Hv_Item Item, Hv_Region region)
{
  int		i;
  char		text[40];
  Hv_String	*Hvtext;
  
  Hv_DrawLine(xc-14*nsize, yc+60+14*nsize, xc+14*nsize, yc+60+14*nsize, Hv_black);
  Hv_DrawLine(xc+60+14*nsize, yc-14*nsize, xc+60+14*nsize, yc+14*nsize, Hv_black);
  for (i=0; i<28; i++) {
    Hv_DrawLine(xc-14*nsize+nsize/2+i*nsize, yc+63+14*nsize, xc-14*nsize+nsize/2+i*nsize, yc+57+14*nsize, Hv_black);
    Hv_DrawLine(xc+63+14*nsize, yc-14*nsize+nsize/2+i*nsize, xc+57+14*nsize, yc-14*nsize+nsize/2+i*nsize, Hv_black);
  }
  Hvtext = Hv_CreateString("0");		   
  Hv_StringDraw(xc+14*nsize-9, yc+14*nsize+65, Hvtext, NULL);
  Hv_StringDraw(xc+14*nsize+68, yc+14*nsize-12, Hvtext, NULL);
  Hvtext = Hv_CreateString("5");		   
  Hv_StringDraw(xc+14*nsize-69, yc+14*nsize+65, Hvtext, NULL);
  Hv_StringDraw(xc+14*nsize+68, yc+14*nsize-72, Hvtext, NULL);
  Hvtext = Hv_CreateString("10");		   
  Hv_StringDraw(xc+14*nsize-132, yc+14*nsize+65, Hvtext, NULL);
  Hv_StringDraw(xc+14*nsize+68, yc+14*nsize-132, Hvtext, NULL);
  Hvtext = Hv_CreateString("15");		   
  Hv_StringDraw(xc+14*nsize-192, yc+14*nsize+65, Hvtext, NULL);
  Hv_StringDraw(xc+14*nsize+68, yc+14*nsize-192, Hvtext, NULL);
  Hvtext = Hv_CreateString("20");		   
  Hv_StringDraw(xc+14*nsize-253, yc+14*nsize+65, Hvtext, NULL);
  Hv_StringDraw(xc+14*nsize+68, yc+14*nsize-252, Hvtext, NULL);
  Hvtext = Hv_CreateString("25");		   
  Hv_StringDraw(xc+14*nsize-313, yc+14*nsize+65, Hvtext, NULL);
  Hv_StringDraw(xc+14*nsize+68, yc+14*nsize-312, Hvtext, NULL);
  Hv_DestroyString(Hvtext);
}

/* - Connects photons in LGB to individual photon data window - */

void DrawConnect(Hv_Item Item, Hv_Region region)
{
  int		h, i, div;

  if (nhit > 5) div=nhit;
  else div=5;
  
  for (i=0; i<nhit; i++) {
    h=(Box4->area->height/div)*((float)i+0.5);
    Hv_DrawLine(Box4->area->left, Box4->area->top+h, xc-(int)(x[i]*ppc), yc-(int)(y[i]*ppc), Hv_black);
  }
}

/* - draw window with data (position & energy) for each photon
     hitting the detector - */

void Photon(Hv_Item Item, Hv_Region region)
{
  Hv_Rect	rect, area;
  int		h, i, j, div, n, length, shift=0;
  char		text[50];
  Hv_String	*Hvtext;

  Hv_CopyRect(&area, Hv_GetItemArea(Item));
  Hv_SetRect(&rect, area.left+3, area.top+3, area.width-6, area.height-6);
  Hv_FillRect(&rect, Hv_gray12);
  
  if (nhit > 5) {
    div=nhit;
    shift=(div-5)*1.8;
  } else
    div=5;
  
  /* - loop over the number of photons - */
  for (i=0; i < nhit; i++) {
    h=(area.height/div)*i;
    sprintf(text, "%d", i+1);
    Hvtext = Hv_CreateString(text);		   
    Hvtext->font=Hv_timesBold17;   
    Hv_StringDraw(area.left+14, area.top+16+h-shift, Hvtext, NULL);
    sprintf(text, "Energy:  %4.3f GeV    Position X: %3.1f", e[i], x[i]);
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+40, area.top+7+h-shift, Hvtext, NULL);
    sprintf(text, "                               Y: %3.1f", y[i]);		   
    Hvtext = Hv_CreateString(text);		   
    Hv_StringDraw(area.left+40, area.top+23+h-2*shift, Hvtext, NULL);
    Hv_DestroyString(Hvtext);

    /* - draws the  energy bars - */
    length=(int)(e[i]*50.0);
    Hv_SetRect(&rect, area.left+40, area.top+26+h-2*shift, length, 9);
    Hv_FillRect(&rect, Hv_red);
    Hv_DrawLine(area.left+40, area.top+30+h-2*shift, area.left+220, area.top+30+h-2*shift, Hv_black);
    Hv_DrawLine(area.left+215, area.top+27+h-2*shift, area.left+220, area.top+30+h-2*shift, Hv_black);
    Hv_DrawLine(area.left+215, area.top+33+h-2*shift, area.left+220, area.top+30+h-2*shift, Hv_black);
    for (j=0; j<7; j++)
      Hv_DrawLine(area.left+40+j*25, area.top+25+h-2*shift, area.left+40+j*25, area.top+35+h-2*shift, Hv_black);
    if (i != 0)
      Hv_DrawLine(area.left+3, area.top+h, area.right-3, area.top+h, Hv_black);
  }
  h=(area.height/div)*i;
  Hv_DrawLine(area.left+3, area.top+h, area.right-3, area.top+h, Hv_black);
}

/* - draws window with data for the LGBs that has been marked by the mouse - */

void Marked(Hv_Item Item, Hv_Region region)
{
  Hv_Rect	rect, area;
  int		i, j, length;
  char		text[50];
  Hv_String	*Hvtext;

  Hv_CopyRect(&area, Hv_GetItemArea(Item));
  Hv_SetRect(&rect, area.left+3, area.top+3, area.width-6, area.height-6);
  Hv_FillRect(&rect, Hv_gray12);
    
  sprintf(text, "Kinematics of %1d marked photons:", nmarked);
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left+5, area.top+5, Hvtext, NULL);
  sprintf(text, "Total energy:     %4.4f GeV", etphot);
  Hv_DrawLine(area.left, area.top+24, area.right, area.top+24, Hv_black);
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left+5, area.top+30, Hvtext, NULL);
  sprintf(text, "Invariant mass:   %4.4f GeV", miphot);
  Hvtext = Hv_CreateString(text);		   
  Hv_StringDraw(area.left+5, area.top+50, Hvtext, NULL);
  sprintf(text, "t-Invariant:     %4.4f GeV\\s\\-\\^2", tiphot);
  Hvtext = Hv_CreateString(text);		   
  Hv_DrawCompoundString(area.left+5, area.top+70, Hvtext);
  Hv_DestroyString(Hvtext);
    
  /* - draws bars to display data - 
  length=etphot/100;
  Hv_SetRect(&rect, area.right-140, area.top+27, length, 9);
  Hv_FillRect(&rect, Hv_red);
  Hv_DrawLine(area.right-140, area.top+31, area.right-10, area.top+31, Hv_black);
  Hv_DrawLine(area.right-15, area.top+28, area.right-10, area.top+31, Hv_black);
  Hv_DrawLine(area.right-15, area.top+34, area.right-10, area.top+31, Hv_black);
  for (j=0; j<6; j++)
    Hv_DrawLine(area.right-140+j*20, area.top+25, area.right-140+j*20, area.top+37, Hv_black);
 */
}

/* - markes/unmarkes photons clicked upon with the mouse - */

void MarkPhoton(Hv_Event hvevent)
{
  int		i, xs, ys, rad;
  float		dist, esum=0.0, xsum=0.0, ysum=0.0, zsum=0.0;

  for (i=0; i < nhit; i++) {
    dist=(worldx-x[i])*(worldx-x[i])+(worldy-y[i])*(worldy-y[i]);
    /* - if the mouse was inside the circle... (compared in cm) - */
    if (dist < (e[i]*12.0/ppc)*(e[i]*12.0/ppc)) {
      if (mark[i] == 0) mark[i]=1;
      else if (mark[i] == 1) mark[i]=0;
      Hv_DrawItem(LeadGlassItem, NULL);
    }
  }
  nmarked=0;
  etphot=0.0;
  miphot=0.0;
  tiphot=0.0;
  for (i=0; i < nhit; i++) {
    if (mark[i] != 0) {
      ++nmarked;
      etphot += e[i];
      esum += e[i];
      xsum += px[i];
      ysum += py[i];
      zsum += pz[i];      
    }
  }
  if (nmarked > 1)
    miphot = sqrt(esum*esum-xsum*xsum-ysum*ysum-zsum*zsum);
  esum -= ebeam;
  zsum -= ebeam;
  tiphot = (esum*esum-xsum*xsum-ysum*ysum-zsum*zsum);
      
  Hv_DrawItem(Box5, NULL);
}
