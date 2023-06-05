/* The event display here represented was designed, coded & documented by   */
/* Troels Petersen under the supervision of Elton Smith, with help from     */
/*   Dave  Heddle and James Muguira.					    */
/* Comments can be found around more complicated code leaving the basics    */
/* uncommented. Upon problems and/or questions, please contact Troels       */
/* Petersen at: <petersen@cebaf.gov> or <petersen@fys.ku.dk>		    */

#include "Event.h"
#include <Hv.h>

static void WindowTile();
static void DrawTile();
static void CenterWrite();
static void OptionsControl();
static void PedestalsOff();
static void PedestalsOn();
static void MAMOn();
static void MAMOff();

Hv_Widget OptionsMenu, optnoped, optped, optmam, optnomam, optmamthresh, optnext, optfile;     

void Init(void)
{
  Hv_View	view;
  int		left, top;
  
  WindowTile();
  
  geo[0]=5;
  geo[1]=6;
  geo[2]=8;
  geo[3]=9;
  geo[4]=10;
  geo[5]=11;
  geo[6]=12;
  geo[7]=12;
  geo[8]=13;
  geo[9]=14;
  geo[10]=14;
  geo[11]=14;
  geo[12]=14;
  geo[13]=14;
  geo[14]=14;

  left = 10;
  top  = 68;


  Hv_VaCreateView(&view,
		  Hv_TAG,             TAGGER_VIEW,
		  Hv_LEFT,            125,
		  Hv_TOP,             665,
		  Hv_HOTRECTLMARGIN,  5,
		  Hv_HOTRECTRMARGIN,  5,
		  Hv_HOTRECTBMARGIN,  80,
                  Hv_HOTRECTCOLOR,    Hv_honeydew,
		  Hv_HOTRECTHEIGHT,   150,
		  Hv_HOTRECTWIDTH,    500,
		  Hv_XMAX,            0.4,
		  Hv_YMAX,            0.66,
		  Hv_INITIALIZE,      SetupTagView,
		  Hv_TITLE,           "Photon Tagger",
		  Hv_POPUPCONTROL,    Hv_STANDARDPOPUPCONTROL - Hv_QUICKZOOM,
		  Hv_DRAWCONTROL,     Hv_STANDARDDRAWCONTROL - Hv_ZOOMONRESIZE,
		  Hv_FEEDBACK,        FeedbackTagger,
		  Hv_OPENATCREATION,  True,
		  NULL);

  Hv_VaCreateView(&view,
		  Hv_TAG,             UPV_VIEW,
		  Hv_LEFT,            135,
		  Hv_TOP,             35,
		  Hv_HOTRECTLMARGIN,  15,
		  Hv_HOTRECTRMARGIN,  15,
		  Hv_HOTRECTBMARGIN,  150,
                  Hv_HOTRECTCOLOR,    Hv_black,
		  Hv_HOTRECTHEIGHT,   250,
		  Hv_HOTRECTWIDTH,    200,
		  Hv_XMAX,            0.4,
		  Hv_YMAX,            0.66,
		  Hv_INITIALIZE,      SetupUPVView,  
		  Hv_TITLE,           "Upstream Particle Veto  (UPV)",
		  Hv_POPUPCONTROL,    Hv_STANDARDPOPUPCONTROL - Hv_QUICKZOOM,
		  Hv_DRAWCONTROL,     Hv_STANDARDDRAWCONTROL - Hv_ZOOMONRESIZE,
		  Hv_FEEDBACK,        FeedbackUPV,
		  Hv_OPENATCREATION,  True,
		  NULL);

  Hv_VaCreateView(&view,
		  Hv_TAG,             CPV_VIEW,
		  Hv_LEFT,            75,
		  Hv_TOP,             51,
		  Hv_HOTRECTLMARGIN,  15,
		  Hv_HOTRECTRMARGIN,  15,
		  Hv_HOTRECTBMARGIN,  130,
                  Hv_HOTRECTCOLOR,    Hv_black,
		  Hv_HOTRECTHEIGHT,   430,
		  Hv_HOTRECTWIDTH,    400,
		  Hv_XMAX,            0.4,
		  Hv_YMAX,            0.66,
		  Hv_INITIALIZE,      SetupCPVView,  
		  Hv_TITLE,           "Charged Particle Veto  (CPV)",
		  Hv_POPUPCONTROL,    Hv_STANDARDPOPUPCONTROL - Hv_QUICKZOOM,
		  Hv_DRAWCONTROL,     Hv_STANDARDDRAWCONTROL - Hv_ZOOMONRESIZE,
		  Hv_FEEDBACK,        FeedbackCPV,
		  Hv_OPENATCREATION,  True,
		  NULL);

  Hv_VaCreateView(&view,
		  Hv_TAG,             BSD_VIEW,
		  Hv_LEFT,            476,
		  Hv_TOP,             51,
		  Hv_HOTRECTLMARGIN,  15,
		  Hv_HOTRECTRMARGIN,  15,
		  Hv_HOTRECTBMARGIN,  130,
                  Hv_HOTRECTCOLOR,    Hv_black,
		  Hv_HOTRECTHEIGHT,   250,
		  Hv_HOTRECTWIDTH,    576,
		  Hv_XMAX,            0.4,
		  Hv_YMAX,            0.66,
       		  Hv_INITIALIZE,      SetupBSDView,  
		  Hv_TITLE,           "Barrel Scintillator Detector (BSD)",
		  Hv_POPUPCONTROL,    Hv_STANDARDPOPUPCONTROL - Hv_QUICKZOOM,
		  Hv_DRAWCONTROL,     Hv_STANDARDDRAWCONTROL - Hv_ZOOMONRESIZE,
		  Hv_FEEDBACK,        FeedbackBSD,
		  Hv_OPENATCREATION,  True,
		  NULL);
  
  Hv_VaCreateView(&view,
		  Hv_TAG,             ED_VIEW,
		  Hv_LEFT,            left,
		  Hv_TOP,             top,
		  Hv_HOTRECTHMARGIN,  85,
		  Hv_HOTRECTWIDTH,    510,
		  Hv_HOTRECTHEIGHT,   510,
                  Hv_HOTRECTCOLOR,    Hv_gray12,
		  Hv_XMAX,            0.4,
		  Hv_YMAX,            0.66,
		  Hv_WIDTH,           900,
		  Hv_HEIGHT,          580,
		  Hv_INITIALIZE,      SetupEdView,
		  Hv_TITLE,           "Rad Phi",
		  Hv_POPUPCONTROL,    Hv_STANDARDPOPUPCONTROL - Hv_QUICKZOOM,
		  Hv_DRAWCONTROL,     Hv_STANDARDDRAWCONTROL - Hv_ZOOMONRESIZE,
		  Hv_FEEDBACK,        Feedback,
		  Hv_OPENATCREATION,  True,
		  NULL);
		  
  Hv_AddFunctionKeyCallback(1, ReadForward);
  
  Hv_VaCreateMenu(&OptionsMenu, Hv_TITLE, "Options" , NULL);
  
  optnoped = Hv_VaCreateMenuItem(OptionsMenu,
		       	      Hv_LABEL,       "pedestals off" ,
			      Hv_CALLBACK,     PedestalsOff,
			      Hv_ACCELERATORCHAR, 'O',
			      NULL);
  
  optped = Hv_VaCreateMenuItem(OptionsMenu,
			      Hv_LABEL,       "pedestals on" ,
			      Hv_CALLBACK,     PedestalsOn,
			      Hv_ACCELERATORCHAR, 'S',
                              NULL);

  optmam =  Hv_VaCreateMenuItem(OptionsMenu,
			      Hv_LABEL,       "display MAM Energy" ,
			      Hv_CALLBACK,     MAMOn,
			      Hv_ACCELERATORCHAR, 'M',
                              NULL);

  optnomam = Hv_VaCreateMenuItem(OptionsMenu,
			      Hv_LABEL,       "no MAM display" ,
			      Hv_CALLBACK,     MAMOff,
			      Hv_ACCELERATORCHAR, 'C',
                              NULL);
 
  optmamthresh = Hv_VaCreateMenuItem(OptionsMenu,
			      Hv_LABEL,       "Minimum MAM Energy" ,
			      Hv_CALLBACK,     MAMSelect,
			      Hv_ACCELERATORCHAR, 'E',
                              NULL);
 
  optnext = Hv_VaCreateMenuItem(OptionsMenu,
			      Hv_LABEL,       "Next Event" ,
			      Hv_CALLBACK,     ReadForward,
			      Hv_ACCELERATORCHAR, 'N',
                              NULL);

  optfile = Hv_VaCreateMenuItem(OptionsMenu,
			      Hv_LABEL,       "New Event File" ,
			      Hv_CALLBACK,     NewFile,
			      Hv_ACCELERATORCHAR, 'F',
                              NULL);
 
  Hv_SetSensitivity(optped, False);   /* default is pedestal subtraction */
  Hv_SetSensitivity(optmam, False);   /* default is display MAM Energy */

}

static void  PedestalsOff(Hv_Widget optnoped)
{
  Hv_SetSensitivity(optnoped, False);
  Hv_SetSensitivity(optped, True);
  pedestals=0;
}

static void  PedestalsOn(Hv_Widget optped)
{
  Hv_SetSensitivity(optped, False);
  Hv_SetSensitivity(optnoped, True);
  pedestals=1;
}

static void  MAMOff(Hv_Widget optnoped)
{
  Hv_SetSensitivity(optnomam, False);
  Hv_SetSensitivity(optmam, True);
  mam_display=0;
}

static void  MAMOn(Hv_Widget optnoped)
{
  Hv_SetSensitivity(optmam, False);
  Hv_SetSensitivity(optnomam, True);
  mam_display=1;
}

static void WindowTile(void)
{
  Hv_tileX = 330;
  Hv_tileY = 80;
  Hv_TileDraw = DrawTile;
}

static void DrawTile()
{
  short      xc, y, color;
  Hv_Rect    rect;
  Hv_Region  region;

  Hv_SetRect(&rect, 0, 0, (short)Hv_tileX, (short)Hv_tileY);
  region = Hv_RectRegion(&rect);
  Hv_SetClippingRegion(region);

  Hv_FillRect(&rect, Hv_canvasColor);
  Hv_DrawLine(rect.left,  rect.bottom, rect.left,  rect.top,    Hv_gray11);
  Hv_DrawLine(rect.left,  rect.top,    rect.right, rect.top,    Hv_gray11);
  Hv_DrawLine(rect.right-1, rect.top,    rect.right-1, rect.bottom-1, Hv_gray5);
  Hv_DrawLine(rect.right-1, rect.bottom-1, rect.left,  rect.bottom-1, Hv_gray5);

  xc = ((short)Hv_tileX)/2;
  y = 14;

  if (Hv_canvasColor == Hv_blue)
    color = Hv_wheat;
  else
    color = Hv_blue;

  CenterWrite("Troels C. Petersen", xc, &y, Hv_helvetica14, color);
  CenterWrite("TJNAF", xc, &y, Hv_helvetica14, color);
  CenterWrite("Radiative Phi Decay", xc, &y, Hv_helvetica14, color);

  Hv_DestroyRegion(region);
}

static void CenterWrite(str, xc, y, font, color)

char           *str;
short          xc;
short          *y;
short          font;
short          color;

{
  Hv_String      *hvstr;
  short          w, h;

  hvstr = Hv_CreateString(str);

  hvstr->font = font;
  hvstr->strcol = color; 

  Hv_CompoundStringDimension(hvstr, &w, &h);
  Hv_DrawCompoundString(xc - w/2, *y, hvstr);

  *y += h;
  Hv_DestroyString(hvstr);
}
