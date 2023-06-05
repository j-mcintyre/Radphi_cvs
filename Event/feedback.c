/* The event display here represented was designed, coded & documented by   */
/* Troels Petersen under the supervision of Elton Smith, with help from     */
/*   Dave  Heddle and James Muguira.					    */
/* Comments can be found around more complicated code leaving the basics    */
/* uncommented. Upon problems and/or questions, please contact Troels       */
/* Petersen at: <petersen@cebaf.gov> or <petersen@fys.ku.dk>		    */

#include  "Event.h"
#include <triggerType.h>

int	lgx, lgy, lgn;

void Feedback(Hv_View View, Hv_Point pp)
{
  float	dist, theta;
  char	text[50];
  int	posx, posy, dx, dy;

  worldx = (float)(xc-pp.x)/(float)nsize*4.0;
  worldy = (float)(yc-pp.y)/(float)nsize*4.0;
  dist = sqrt(worldx*worldx+worldy*worldy);
  theta = atan2(dist,100.0);

  sprintf(text, "Real World X: %.2f cm", worldx);
  Hv_ChangeFeedback(View, WORLDX, text);

  sprintf(text, "Real World Y: %.2f cm", worldy);
  Hv_ChangeFeedback(View, WORLDY, text);

  sprintf(text, "Distance from beam: %.2f cm", dist);
  Hv_ChangeFeedback(View, DIST, text);

  sprintf(text, "Angle from beam: %.3f rad", theta);
  Hv_ChangeFeedback(View, THETA, text);

  /* - geometry conversion (from symmetry) for the IF-statement - */
  posx=(pp.x-xc)/nsize;
  if (posx < 0) posx=-posx;
  dy=geo[13-posx]*nsize;
  posy=(pp.y-yc)/nsize;
  if (posy < 0) posy=-posy;
  dx=geo[13-posy]*nsize;
  
  /* - if the mouse is at a LGB, then... - */
  if ((xc-dx<pp.x) && (pp.x<xc+dx) && (yc-dy<pp.y) && (pp.y<yc+dy)) { /* in the LGD */

    lgx = 27-(pp.x-xp)/nsize;
    lgy = 27-(pp.y-yp)/nsize;
    lgn = lgx + lgy*28;
    if (((lgx==13) || (lgx==14)) && ((lgy==13) || (lgy==14))) {
      lgx=0; lgy=0; lgn=0;
    }  
    sprintf(text, "ID number: %4d   Lead Glass X: %3d (column)", lgn, lgx);
    Hv_ChangeFeedback(View, POSS, text);

    if (calibrated && trigbit==TRIG_DATA) {
      sprintf(text, "Energy (MeV):%4d Lead Glass Y: %3d (row)", (int)lg[lgx][lgy], lgy);
      Hv_ChangeFeedback(View, INFO, text);
    }
    else {
       sprintf(text, "ADC:    %4d      Lead Glass Y: %3d (row)", (int)lg[lgx][lgy], lgy);
       Hv_ChangeFeedback(View, INFO, text);
    }
  } else {
    Hv_ChangeFeedback(View, POSS, "                                   ");
    Hv_ChangeFeedback(View, INFO, "                                   ");
  }  
}

void FeedbackTagger(Hv_View View, Hv_Point pp)
{

  int   dx, dy, id;

  /*  average energy of photon w.r.t. original beam energy */
  float  tag_energy[19] = {0.9490, 0.93975, 0.9288, 0.9183,
                           0.9096, 0.9006, 0.89195, 0.8836,
                           0.87445, 0.8648, 0.85475, 0.8441,
                           0.8329, 0.8219, 0.8129, 0.80565,
                           0.79795, 0.7897, 0.78105};
 

  worldx = (float)(pp.x)/(float)nsize*4.0;
  worldy = (float)(pp.y)/(float)nsize*4.0;

  dx = pp.x - xctag ;
  dy = pp.y - yctag ; 

  id = dx/nsize;    /* note: hardware numbering starts at 1 */

  if (dy >= 0 && dy < nsize*2 && id < 19 ) {
    /* continue */
  }
  else if (dy <0 && dy > -nsize*2 && id < 19 ) {
    id += 19;
  }
  else {
    id = -1;
  }

  if (id >= 0) {
    char ctext[50];
    char etext[50];
    char ttext[50];
    float photon_energy;
    if (id < 19) {
      sprintf(ctext, "Tagger channel:  %.1i Left",id+1);
      photon_energy = beamEnergy * tag_energy[id];
    }
    else {
      sprintf(ctext, "Tagger channel:  %.1i Right",id-18);
      photon_energy = beamEnergy * tag_energy[id-19];
    }
    sprintf(etext, "Photon Energy: %.2f GeV", photon_energy);
    Hv_ChangeFeedback(View,8,ctext);
    Hv_ChangeFeedback(View,9,"                                   ");
    Hv_ChangeFeedback(View,10,"                                   ");
    Hv_ChangeFeedback(View,22,etext);
    if (calibrated) {
      int i;
      int line = 9;
      for (i = 0; i < tagger_times[id].multicount; i++) {
        float time = tagger_times[id].time[i];
        if (i == 0) {
          sprintf(ttext, "Time list:    %6.1f ns", time);
        }
        else {
          sprintf(ttext, "              %6.1f ns", time);
        }
        Hv_ChangeFeedback(View,line++,ttext);
      }
    }
    else {
      int i;
      int line = 9;
      for (i = 0; i < tagger_times[id].multicount; i++) {
        int tdc = tagger_times[id].time[i];
        if (i == 0) {
          sprintf(ttext, "TDC list:    %4i", tdc);
        }
        else {
          sprintf(ttext, "             %4i", tdc);
        }
        Hv_ChangeFeedback(View,line++,ttext);
      }
    }
  }
  else {
    Hv_ChangeFeedback(View,  8, "                                   ");
    Hv_ChangeFeedback(View,  9, "                                   ");
    Hv_ChangeFeedback(View, 10, "                                   ");
    Hv_ChangeFeedback(View, 22, "                                   ");
  }
}

void FeedbackCPV(Hv_View View, Hv_Point pp)
{  
  float ppi = 7.9558;
  int 	id;
  worldx = pp.x - xctcpv;
  worldy = pp.y - yctcpv;
 
  if (((worldy > 15.875*ppi) && (worldy < 22.625*ppi)) &&
       (fabs(worldx) < 17.188*ppi)) { 
    id = (worldx > 0)? 0 : 1;
  }
  else if (((worldy > 12.375*ppi) && (worldy < 15.875*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) {
    id = (worldx > 0)? 2 : 3;
  }
  else if (((worldy > 8.875*ppi) && (worldy < 12.375*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) { 
    id = (worldx > 0)? 4 : 5;
  }
  else if (((worldy > 5.375*ppi) && (worldy < 8.875*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) {
    id = (worldx > 0)? 6 : 7;
  }
  else if (((worldy > 3.125*ppi) && (worldy < 5.375*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) { 
    id = (worldx > 0)? 8 : 9;
  }
  else if (((worldy > 1.875*ppi) && (worldy < 3.125*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) { 
    id = (worldx > 0)? 10 : 11;
  }
  else if (((worldy > 0.625*ppi) && (worldy < 1.875*ppi)) &&
           ((fabs(worldx) > 1.626*ppi) && (fabs(worldx) < 22.188*ppi))) {
    id = (worldx > 0)? 12 : 13;
  }
  else if (((worldy < -15.875*ppi) && (worldy > -22.625*ppi)) &&
            (fabs(worldx) < 17.188*ppi)) { 
    id = (worldx > 0)? 28 : 29;
  } 
  else if (((worldy < -12.375*ppi) && (worldy > -15.875*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) { 
    id = (worldx > 0)? 26 : 27;
  }
  else if (((worldy < -8.875*ppi) && (worldy > -12.375*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) { 
    id = (worldx > 0)? 24 : 25;
  }
  else if (((worldy < -5.375*ppi) && (worldy > -8.875*ppi)) &&
            (fabs(worldx) < 22.188*ppi)){ 
    id = (worldx > 0)? 22 : 23;
  } 
  else if (((worldy < -3.125*ppi) && (worldy > -5.375*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) { 
    id = (worldx > 0)? 20 : 21;
  }
  else if (((worldy < -1.875*ppi) && (worldy > -3.125*ppi)) &&
            (fabs(worldx) < 22.188*ppi)) { 
    id = (worldx > 0)? 18 : 19;
  }
  else if (((worldy < -0.625*ppi) && (worldy > -1.875*ppi)) &&
           ((fabs(worldx) < 22.188*ppi) && (fabs(worldx) > 1.626*ppi))) { 
    id = (worldx > 0)? 16 : 17;
  }
  else if (((worldy < .625*ppi) && (worldy > -.625*ppi)) &&
           ((fabs(worldx) > 1.626*ppi) && (fabs(worldx) < 22.188*ppi))) {
    id = (worldx > 0)? 14 : 15;
  }
  else{
    id = -1;
  }

  if (id >= 0) {
    char ctext[50];
    char etext[50];
    char ttext[50];
    sprintf(ctext, "CPV channel: %i", id);
    Hv_ChangeFeedback(View,10,ctext);
    Hv_ChangeFeedback(View,11,"                       ");
    Hv_ChangeFeedback(View,12,"                       ");
    Hv_ChangeFeedback(View,13,"                       ");
    Hv_ChangeFeedback(View,14,"                       ");
    if (calibrated) {
      float energy = cpv[id];
      if (energy > 0) {
        sprintf(etext, "deltaE value: %i KeV", (int)energy);
        Hv_ChangeFeedback(View,11,etext);
      }
      if (cpv_times[id].multicount > 0) {
        int i;
        int line = 12;
	for (i = 0;i < cpv_times[id].multicount; i++) {
	  float time = (float)cpv_times[id].time[i];
          if (i == 0) {
	    sprintf(ttext, "Time list: %6.1f ns",time);
          }
          else {
	    sprintf(ttext, "           %6.1f ns",time);
          }
          Hv_ChangeFeedback(View,line++,ttext);
	}
      }
    }
    else {
      int adc = cpv[id];
      if (adc > 0) {
        sprintf(etext, "ADC value  : %4i", adc);
        Hv_ChangeFeedback(View,11,etext);
      }
      if (cpv_times[id].multicount > 0){
        int i;
        int line = 12;
	for (i = 0;i < cpv_times[id].multicount; i++) {
	  int tdc = (float)cpv_times[id].time[i];
          if (i == 0) {
	    sprintf(ttext, "TDC list:  %4i",tdc);
          }
          else {
	    sprintf(ttext, "           %4i",tdc);
          }
          Hv_ChangeFeedback(View,line++,ttext);
        }
      }
    }
  }
  else {
    Hv_ChangeFeedback(View,10,"               ");
    Hv_ChangeFeedback(View,11,"               ");
    Hv_ChangeFeedback(View,12,"               ");
    Hv_ChangeFeedback(View,13,"               ");
    Hv_ChangeFeedback(View,14,"               ");
  }
}

void FeedbackBSD(Hv_View View, Hv_Point pp) 
{
  int 	id_s, id_r, id_l;
  int	a,b,c,d,rightx,righty,leftx,lefty;
  float time_s, time_r, time_l;
  int	coordx, coordy;
  int 	SW = 24; /* distance between straight counters in pixels */
  int 	LW = 48; /* distance between left or right counters in pixels */

  /* obtain the world coordinates from the mouse location */

  Hv_LocalToWorld(View, &worldx, &worldy, pp.x, pp.y);

  worldx = (pp.x - xctbsd);
  worldy = (pp.y - yctbsd);

  coordx = worldx + 288;
  coordy = worldy - 96;

  /* find straight counter where mouse is located */

  if ((worldx > -288 && worldx < 288) &&
      (worldy >  -96 && worldy <  96)) {
    id_s = (coordx - coordx%SW) / SW;

    /* do shear transformation to find right counter where mouse is located */
    /* transformation = [(a,c),(b,d)] */
    
    a=1;
    b=1;
    c=0;
    d=1;
    
    rightx = a*coordx + b*coordy;
    righty = c*coordx + d*coordy;

    if (rightx > 0) {
      id_r = (rightx - rightx%LW) / LW;
    }
    else {
      rightx = abs(rightx);
      id_r = -((rightx - rightx%LW) / LW) + 11;
    }

    /* do shear transformation to find left counter where mouse is located */
    /* transformation = [(a,c),(b,d)] */
    
    a=1;
    b=-1;
    c=0;
    d=1;
    
    leftx = a*coordx + b*coordy;
    lefty = c*coordx + d*coordy;
    
    if (leftx < 576) {
      id_l = (leftx - leftx%LW) / LW;
    }
    else {
      leftx = abs(leftx);
      id_l = ((leftx - leftx%LW) / LW) - 12;
    }
  }
  else {
    id_s = -1;
    id_r = -1;
    id_l = -1;
  }

  if (id_s >= 0) {
    char setext[50];
    char retext[50];
    char letext[50];
    char sttext[50];
    char rttext[50];
    char lttext[50];
    sprintf(setext,"BSD Counter: S%2d  deltaE:       ",id_s);
    sprintf(retext,"           : R%2d        :       ",id_r);
    sprintf(letext,"           : L%2d        :       ",id_l);
    sprintf(sttext,"BSD Counter: S%2d    time:        ",id_s);
    sprintf(rttext,"             R%2d        :        ",id_r);
    sprintf(lttext,"             L%2d        :        ",id_l);
    if (calibrated) {
      int i;
      float se = bsd_s[id_s];
      float sr = bsd_r[id_r];
      float sl = bsd_l[id_l];
      if (se > 0) {
        sprintf(setext,"BSD Counter: S%2d  deltaE: %i KeV",id_s,(int)se);
      }
      if (sr > 0) {
        sprintf(retext,"           : R%2d        : %i KeV",id_r,(int)sr);
      }
      if (sl > 0) {
        sprintf(letext,"           : L%2d        : %i KeV",id_l,(int)sl);
      }
      for (i = 0; i < bsd_s_times[id_s].multicount; i++) {
        float st = bsd_s_times[id_s].time[i];
        sprintf(sttext,"BSD Counter: S%2d    time: %.1f ns",id_s,st);
      }
      for (i = 0; i < bsd_r_times[id_r].multicount; i++) {
        float rt = bsd_r_times[id_r].time[i];
        sprintf(rttext,"             R%2d        : %.1f ns",id_r,rt);
      }
      for (i = 0; i < bsd_l_times[id_l].multicount; i++) {
        float lt = bsd_l_times[id_l].time[i];
        sprintf(lttext,"             L%2d        : %.1f ns",id_l,lt);
      }
    }
    else {
      int i;
      int sadc = bsd_s[id_s];
      int radc = bsd_r[id_r];
      int ladc = bsd_l[id_l];
      if (sadc > 0) {
        sprintf(setext,"BSD Counter: S%2d     ADC: %i",id_s,sadc);
      }
      if (radc > 0) {
        sprintf(retext,"           : R%2d        : %i",id_r,radc);
      }
      if (ladc > 0) {
        sprintf(letext,"           : L%2d        : %i",id_l,ladc);
      }
      for (i = 0; i < bsd_s_times[id_s].multicount; i++) {
        int stdc = bsd_s_times[id_s].time[i];
        sprintf(sttext,"BSD Counter: S%2d     TDC: %.1f ns",id_s,stdc);
      }
      for (i = 0; i < bsd_r_times[id_r].multicount; i++) {
        int rtdc = bsd_r_times[id_r].time[i];
        sprintf(rttext,"             R%2d        : %.1f ns",id_r,rtdc);
      }
      for (i = 0; i < bsd_l_times[id_l].multicount; i++) {
        int ltdc = bsd_l_times[id_l].time[i];
        sprintf(lttext,"             L%2d        : %.1f ns",id_l,ltdc);
      }
    }
    Hv_ChangeFeedback(View,16,setext);
    Hv_ChangeFeedback(View,17,retext);
    Hv_ChangeFeedback(View,18,letext);
    Hv_ChangeFeedback(View,19,sttext);
    Hv_ChangeFeedback(View,20,rttext);
    Hv_ChangeFeedback(View,21,lttext);

  } /* end if */
  else {
    Hv_ChangeFeedback(View,16,"                               ");
    Hv_ChangeFeedback(View,17,"                               ");
    Hv_ChangeFeedback(View,18,"                               ");
    Hv_ChangeFeedback(View,19,"                     ");
    Hv_ChangeFeedback(View,20,"                     ");
    Hv_ChangeFeedback(View,21,"                     ");
  } /* end else */
  
} /* end FeedbackBSD */

void FeedbackUPV(Hv_View View, Hv_Point pp)
{
  float ppi = 5.54;
  int  	id;
 
  worldx = (pp.x - xctupv);
  worldy = (pp.y - yctupv); 

  if (((worldx > 2.25*ppi) && (worldx < 16.75*ppi)) &&
      ((worldy > 3.37*ppi) && (worldy < 10.47*ppi))) {
    id = 0;
  }
  else if (((worldx > 2.25*ppi) && (worldx < 16.75*ppi)) &&
           ((worldy > -3.37*ppi) && (worldy < 3.37 *ppi))) {
    id = 2;
  }
  else if (((worldx > 2.25*ppi) && (worldx < 16.75*ppi)) &&
           ((worldy < -3.37*ppi) && (worldy > -10.47*ppi))) {
    id = 4;
  }
  else if (((worldx < -2.25*ppi) && (worldx > -16.75*ppi)) &&
           (( worldy > 3.37*ppi) && (worldy < 10.47*ppi))) {
    id = 1;
  }
  else if (((worldx < -2.25*ppi) && (worldx > -16.75*ppi)) &&
           ((worldy > -3.37*ppi) && (worldy < 3.37 *ppi))) {
    id = 3;
  }
  else if (((worldx < -2.25 *ppi) && (worldx > -16.75*ppi)) &&
           ((worldy < -3.37 *ppi) && (worldy > -10.47*ppi))) {
    id = 5;
  }
  else if ( (((worldx > -3.55*ppi) && (worldx < 3.35*ppi)) &&
             ((worldy < -10.65*ppi) && (worldy > -16.75*ppi))) ||
	    (((worldx > -2.25*ppi) && (worldx < 2.25*ppi)) &&
             ((worldy < -2.25*ppi) && (worldy > -10.65*ppi))) ) {
    id = 7;
  }
  else if ( (((worldx > -3.55*ppi) && (worldx < 3.35*ppi)) &&
             ((worldy > 10.65*ppi) && (worldy < 16.75*ppi))) ||
	    (((worldx > -2.25*ppi) && (worldx < 2.25*ppi)) &&
             ((worldy > 2.25*ppi) && (worldy < 10.65*ppi))) ) {
    id = 6;
  }
  else{
    id = -1;
  }

  if (id >= 0) {
    char ctext[50];
    char etext[50];
    char ttext[50];
    sprintf(ctext, "UPV id   : %i",id);
    Hv_ChangeFeedback(View,13,ctext);
    Hv_ChangeFeedback(View,14,"                               ");
    Hv_ChangeFeedback(View,15,"                               ");
    Hv_ChangeFeedback(View,16,"                               ");
    Hv_ChangeFeedback(View,17,"                               ");
    if (calibrated) {
      float energy = upv[id];
      sprintf(etext, "deltaE value: %i", (int)energy);
      Hv_ChangeFeedback(View,14,etext);
      if(upv_times[id].multicount > 0){
        int i;
        int line = 15;
        for (i = 0; i < upv_times[id].multicount; i++) {
          float time = upv_times[id].time[i];
          if (i == 0) {
	    sprintf(ttext, "Time list:  %6.1f ns",time);
          }
          else {
	    sprintf(ttext, "            %6.1f ns",time);
          }
          Hv_ChangeFeedback(View,line++,ttext);
	}
      }
    }
    else {
      int adc = upv[id];
      sprintf(etext, "ADC value: %i",adc);
      Hv_ChangeFeedback(View,14,etext);
      if (upv_times[id].multicount > 0){
        int i;
        int line = 15;
        for (i = 0;i < upv_times[id].multicount; i++) {
          int tdc = upv_times[id].time[i];
          if (i == 0) {
	    sprintf(ttext, "First TDC: %4i",tdc);
          }
          else {
	    sprintf(ttext, "           %4i",tdc);
          }
          Hv_ChangeFeedback(View,line++,ttext);
	}
      }
    }
  }
  else {
    Hv_ChangeFeedback(View,13,"           ");
    Hv_ChangeFeedback(View,14,"           ");
    Hv_ChangeFeedback(View,15,"           ");
    Hv_ChangeFeedback(View,16,"           ");
    Hv_ChangeFeedback(View,17,"           ");
  }
}
