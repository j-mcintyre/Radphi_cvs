/* The event display here represented was designed, coded & documented by   */
/* Troels Petersen under the supervision of Elton Smith, with help from     */
/*   Dave  Heddle and James Muguira.					    */
/* Comments can be found around more complicated code leaving the basics    */
/* uncommented. Upon problems and/or questions, please contact Troels       */
/* Petersen at: <petersen@cebaf.gov> or <petersen@fys.ku.dk>		    */

#include "Event.h"

Hv_Item  FeedbackBSDItem, FeedbackCPVItem, FeedbackBGVItem,
         FeedbackTaggerItem, FeedbackItem, FeedbackUPVItem;  
Hv_Item	 Scale, Button, Option, Info;

/* - updates the LeadGlassItem after a change in options - */

static void OptionButtonCB(Hv_Event hvevent)
{
  Hv_OptionButtonCallback(hvevent);
  Hv_DrawItem(LeadGlassItem, NULL);
}

/* - definition (position and size) of all the objects
     (controls and information) in the main (Rad Phi) view       - */

static void   ControlItems(Hv_View View)
{
  Hv_Item	Item;
  Hv_Rect	*hotrect = View->hotrect;
			   
/* - creates a camera button for printing - */

  Hv_StandardTools(View, hotrect->left - 63, hotrect->top - 5, Hv_VERTICAL,
		   Hv_DONOTHINGONRESIZE, Hv_CAMERATOOL, &Item, &Item);

/* Vertical energy scale */

  Box1 = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BOXITEM,
		Hv_LEFT,		hotrect->left - 82,
		Hv_TOP,			hotrect->top + 27,
		Hv_WIDTH,		75,
		Hv_HEIGHT,		485,
		NULL);

  Scale = Hv_VaCreateItem(View,
		Hv_PARENT,		Box1,
		Hv_BALLOON,		(char *)"Energy Scale for the LGBs.",
		Hv_TYPE,		Hv_RAINBOWITEM,
		Hv_ORIENTATION,		Hv_VERTICAL, 
		Hv_LEFT,		Box1->area->left + 28,
		Hv_TOP,			Box1->area->top + 50,
		Hv_EDITABLE,		True,
		Hv_DRAWCONTROL,		Hv_REDRAWHOTRECTWHENCHANGED,
		Hv_WIDTH,		15,
		Hv_NUMLABELS,		5,
		Hv_PRECISION,		0,
		Hv_MINVALUE,		0,
		Hv_MAXVALUE,		LGD_RED_HOT,
		Hv_COLORWIDTH,		5,
		Hv_TITLE,		" Energy\\n(in MeV)",   
			  /*		Hv_TITLE,		" Energy\\n(channels)", */
		Hv_MINCOLOR,		Hv_white,
		Hv_MAXCOLOR,		Hv_red,
		NULL);
			    
/* Control buttons */

  Box2 = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BOXITEM,
		Hv_LEFT,		hotrect->right + 15,
		Hv_TOP,			hotrect->top - 3,
		Hv_WIDTH,		200,
		Hv_HEIGHT,		70,
		NULL);


  Button = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button to go to the next event. If at End-Of-File, a window will tell you so.",
		Hv_PARENT,		Box2,
	       	Hv_SINGLECLICK,		ReadForward,
		Hv_LEFT,		Box2->area->left + 102,
		Hv_TOP,			Box2->area->top + 37,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_white,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"Next event",
		NULL);
			 
  Button = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button for jumping to the any event. If chosen event exceeds the file length, a window will tell you so.",
		Hv_PARENT,		Box2,
		Hv_SINGLECLICK,		ReadGoto,
		Hv_LEFT,		Box2->area->left + 4,
		Hv_TOP,			Box2->area->top + 37,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_wheat,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"Goto event",
		NULL);

  Button = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button for loading file.'*.itape' is the default filter.",
		Hv_PARENT,		Box2,
		Hv_SINGLECLICK,		NewFile,
		Hv_LEFT,		Box2->area->left + 102,
		Hv_TOP,			Box2->area->top + 3,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_gray13,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"File",
		NULL);
  Button = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button for getting input data from Dispatcher for online running.",
		Hv_PARENT,		Box2,
			   /*		Hv_SINGLECLICK,		NewInput, */
		Hv_SINGLECLICK,		Online, 
		Hv_LEFT,		Box2->area->left + 4,
		Hv_TOP,			Box2->area->top + 3,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_gray12,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"online",
			   /*		Hv_TEXT,		"stdin", */
		NULL);

/* Option controls */

  Box3 = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BOXITEM,
		Hv_PLACEMENTITEM,	Box2,
		Hv_RELATIVEPLACEMENT,	Hv_POSITIONRIGHT,
		Hv_PLACEMENTGAP,	5,
		Hv_WIDTH,		160,
		Hv_HEIGHT,		70,
		NULL);

  Option = Hv_VaCreateItem(View,
		Hv_PARENT,		Box3,
		Hv_BALLOON,		(char *)"Match up of the photons to mesons ON/OFF.",
		Hv_TYPE,		Hv_OPTIONBUTTONITEM,
		Hv_CONTROLVAR,		(void *)(&matchup),
		Hv_LEFT,		Box3->area->left + 10,
		Hv_TOP,			Box3->area->top + 6,
		Hv_TEXT,		"Match Up",
		Hv_FONT,		Hv_fixed2,
		Hv_ARMCOLOR,		Hv_red,
		Hv_SINGLECLICK,		OptionButtonCB,
		NULL);

  Option = Hv_VaCreateItem(View,
		Hv_PARENT,		Box3,
		Hv_BALLOON,		(char *)"Scale of the LGB's ID numbers ON/OFF.",
		Hv_TYPE,		Hv_OPTIONBUTTONITEM,
		Hv_CONTROLVAR,		(void *)(&lgscale),
		Hv_LEFT,		Box3->area->left + 10,
		Hv_TOP,			Box3->area->top + 26,
		Hv_TEXT,		"LG scale",
		Hv_FONT,		Hv_fixed2,
		Hv_ARMCOLOR,		Hv_red,
		Hv_SINGLECLICK,		OptionButtonCB,
		NULL);

  Option = Hv_VaCreateItem(View,
		Hv_PARENT,		Box3,
		Hv_BALLOON,		(char *)"Connecting lines between photon circles and corresponding data ON/OFF.",
		Hv_TYPE,		Hv_OPTIONBUTTONITEM,
		Hv_CONTROLVAR,		(void *)(&phodata),
		Hv_LEFT,		Box3->area->left + 10,
		Hv_TOP,			Box3->area->top + 46,
		Hv_TEXT,		"Connect",
		Hv_FONT,		Hv_fixed2,
		Hv_ARMCOLOR,		Hv_red,
		Hv_SINGLECLICK,		OptionButtonCB,
		NULL);

  Button = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button to select trigger types. Allows one to choose to displays of a certain trigger type only. Choosing 0 accepts all triggers",
		Hv_PARENT,		Box3,
		Hv_SINGLECLICK,		TrigSelect,   
		Hv_LEFT,		Box3->area->left + 90,
		Hv_TOP,			Box3->area->top + 10,
		Hv_WIDTH,		60,
		Hv_HEIGHT,		40,
		Hv_FILLCOLOR,		Hv_aquaMarine,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"Trigger\\nSelect",
		NULL);


/* FeedbackItem (information on mouse position) */

  FeedbackItem = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_FEEDBACKITEM,
		Hv_BALLOON,		(char *)"Data from mouse position, when at detector.",
		Hv_PLACEMENTITEM,	Box2,
		Hv_RELATIVEPLACEMENT,	Hv_POSITIONBELOW,
		Hv_PLACEMENTGAP,	5,
		Hv_FEEDBACKDATACOLS,	50,
		Hv_FEEDBACKENTRY,	WORLDX, Hv_fixed2, Hv_black, "",
		Hv_FEEDBACKENTRY,	WORLDY, Hv_fixed2, Hv_black, "",
		Hv_FEEDBACKENTRY,	DIST,   Hv_fixed2, Hv_black, "",
		Hv_FEEDBACKENTRY,	THETA,  Hv_fixed2, Hv_black, "",
		Hv_FEEDBACKENTRY,	Hv_SEPARATORTAG,
		Hv_FEEDBACKENTRY,	POSS, Hv_fixed2, Hv_navyBlue, "",
		Hv_FEEDBACKENTRY,	INFO, Hv_fixed2, Hv_navyBlue, "",
		Hv_FEEDBACKMARGIN,	15,
		NULL);

/* Photon data - drawing routine: Photon */

  Box4 = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BOXITEM,
		Hv_BALLOON,		(char *)"Reconstructed data for each photon. Use numbers or 'Connect' for correspondance.",
		Hv_PLACEMENTITEM,	FeedbackItem,
		Hv_RELATIVEPLACEMENT,	Hv_POSITIONBELOW,
		Hv_PLACEMENTGAP,	5,
		Hv_WIDTH,		FeedbackItem->area->width,
		Hv_HEIGHT,		220,
		Hv_USERDRAW,		Photon,
		NULL);

/* Marked data - drawing routine: Marked */

  Box5 = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BOXITEM,
		Hv_BALLOON,		(char *)"Calculations on marked photons. Mark/unmark by clicking inside the photons circle.",
		Hv_PLACEMENTITEM,	Box4,
		Hv_RELATIVEPLACEMENT,	Hv_POSITIONBELOW,
		Hv_PLACEMENTGAP,	5,
		Hv_WIDTH,		FeedbackItem->area->width,
		Hv_HEIGHT,		90,
		Hv_USERDRAW,		Marked,
		NULL);

/* Programmer information */

  Info = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_TEXTITEM,
		Hv_TEXT,		"Troels C. Petersen, NBI/Jefferson Lab, June 1996",
		Hv_FONT,		Hv_times14,
		Hv_PLACEMENTITEM,	Box5,
		Hv_RELATIVEPLACEMENT,	Hv_POSITIONBELOW,
		Hv_PLACEMENTGAP,	5,
		Hv_WIDTH,		FeedbackItem->area->width,
		Hv_HEIGHT,		15,
		NULL);

/* The detector itself */

  LeadGlassItem = Hv_VaCreateItem(View,
		Hv_USERDRAW,		LeadGlass,
		Hv_BALLOON,		(char *)"The LGD and RPD detectors - not drawn to scale!",
		Hv_LEFT,		hotrect->left + 10,
		Hv_TOP,			hotrect->top + 10,
		Hv_WIDTH,		500,
		Hv_HEIGHT,		500,
		Hv_SINGLECLICK,		MarkPhoton,
		NULL);
}

static void   ControlTagItems(Hv_View View)
{
  Hv_Item	Item;
  Hv_Rect	*hotrect = View->hotrect;
			   
/* - creates a camera button for printing - */

  Hv_StandardTools(View, hotrect->left+340, hotrect->bottom + 10, Hv_VERTICAL,
		   Hv_DONOTHINGONRESIZE, Hv_CAMERATOOL, &Item, &Item);

  Button = Hv_VaCreateItem(View,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button to go to the next event.",
		Hv_SINGLECLICK,		ReadForward,
		Hv_LEFT,		hotrect->left + 230,
		Hv_TOP,			hotrect->bottom + 45,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_white,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"Next event",
		NULL);

  Button = Hv_VaCreateItem(View,
		Hv_TYPE,                Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button to select beam energy.",
		Hv_SINGLECLICK,		BeamEnergySelect,
		Hv_LEFT,		hotrect->left + 230,
		Hv_TOP,			hotrect->bottom + 10,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_yellow,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"Beam Energy\\nSelect",
		NULL);

  FeedbackTaggerItem = Hv_VaCreateItem(View,
                Hv_TYPE,                Hv_FEEDBACKITEM,
                Hv_BALLOON,             (char *)"Data from mouse position, when at tagger.",
		Hv_LEFT,                hotrect->left ,
                Hv_TOP,                 hotrect->bottom + 10 ,
                Hv_WIDTH,               1,
                Hv_HEIGHT,              1, 
                Hv_FEEDBACKDATACOLS,    25,
                Hv_FEEDBACKENTRY,       8,  Hv_fixed2, Hv_darkBlue, "", 
                Hv_FEEDBACKENTRY,       9,  Hv_fixed2, Hv_black, "",
                Hv_FEEDBACKENTRY,       10, Hv_fixed2, Hv_black, "",
	       	Hv_FEEDBACKENTRY,       Hv_SEPARATORTAG,
		Hv_FEEDBACKENTRY,       22, Hv_fixed2, Hv_navyBlue, "",
                Hv_FILLCOLOR,           Hv_white,
                NULL);


  TaggerItem = Hv_VaCreateItem(View,
		Hv_USERDRAW,		Tagger,
		Hv_BALLOON,		(char *)"The photon Tagger detectors - schematic. The color represents relative time values",
		Hv_LEFT,		hotrect->left,
		Hv_TOP,			hotrect->top,
		Hv_WIDTH,		500,
		Hv_HEIGHT,		150,
		NULL);
}

static void   ControlCPVItems(Hv_View View)
{
  Hv_Item	Item;
  Hv_Rect	*hotrect = View->hotrect;
			
  Button = Hv_VaCreateItem(View,
		Hv_RESIZEPOLICY,        Hv_SHIFTYONRESIZE,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button to go to the next event.",
		Hv_SINGLECLICK,		ReadForward,
		Hv_LEFT,		hotrect->left + 233,
		Hv_TOP,			hotrect->bottom + 50,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_white,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"Next event",
		NULL);

  Box6 = Hv_VaCreateItem(View,
                Hv_RESIZEPOLICY,        Hv_SHIFTYONRESIZE,
		Hv_TYPE,		Hv_BOXITEM,
		Hv_LEFT,		hotrect->left,
		Hv_TOP,			hotrect->bottom + 20,
		Hv_WIDTH,		213,
		Hv_HEIGHT,		100,
		Hv_FILLCOLOR,		Hv_gray14,
		NULL);
		

 
  FeedbackCPVItem = Hv_VaCreateItem(View,
                Hv_RESIZEPOLICY,        Hv_SHIFTYONRESIZE,
		Hv_TYPE,		Hv_FEEDBACKITEM,
		Hv_BALLOON,		(char *)"Data from mouse position, when at cpv.",
		Hv_PARENT,		Box6,
		Hv_LEFT,		Box6->area->left ,
		Hv_TOP,			Box6->area->top ,
		Hv_WIDTH,               213,
                Hv_HEIGHT,              100, 
                Hv_FEEDBACKDATACOLS,    25,
                Hv_FEEDBACKENTRY,       10, Hv_fixed2, Hv_darkBlue, "",
                Hv_FEEDBACKENTRY,       11, Hv_fixed2, Hv_black, "",
                Hv_FEEDBACKENTRY,       12, Hv_fixed2,  Hv_black, "", 
                Hv_FEEDBACKENTRY,       13, Hv_fixed2,  Hv_black, "", 
                Hv_FEEDBACKENTRY,       14, Hv_fixed2,  Hv_black, "", 
	       	Hv_FEEDBACKMARGIN,      20,
                Hv_FILLCOLOR,           Hv_white,
                NULL);
	    

    CPVItem = Hv_VaCreateItem(View,
		Hv_USERDRAW,		CPV,
		Hv_BALLOON,		(char *)"The CPV data- not drawn to scale! Solid colors represent ADC hits. Vertical stripes represent an ADC hit within the TDC window",
		Hv_LEFT,		hotrect->left,
		Hv_TOP,			hotrect->top,
		Hv_WIDTH,		400,
		Hv_HEIGHT,		430,
		NULL);
}

static void   ControlBSDItems(Hv_View View)
{
  Hv_Item       Item;
  Hv_Rect       *hotrect = View->hotrect;

  Button = Hv_VaCreateItem(View,
                Hv_RESIZEPOLICY,       Hv_SHIFTYONRESIZE,
		Hv_TYPE,               Hv_BUTTONITEM,
	        Hv_BALLOON,            (char *)"Button to go to the next event.",
	        Hv_SINGLECLICK,        ReadForward,
                Hv_LEFT,               hotrect->left + 300,
                Hv_TOP,                hotrect->bottom + 50,
                Hv_WIDTH,              96,
	        Hv_HEIGHT,             30,
                Hv_FILLCOLOR,          Hv_white,
	        Hv_FONT,               Hv_fixed,
	        Hv_TEXT,               "Next event",
                NULL);

  Box6 = Hv_VaCreateItem(View,
		Hv_RESIZEPOLICY,       Hv_SHIFTYONRESIZE,
                Hv_TYPE,               Hv_BOXITEM,
		Hv_LEFT,               hotrect->left,
                Hv_TOP,                hotrect->bottom +20,
                Hv_WIDTH,              283,
                Hv_HEIGHT,             100,
                Hv_FILLCOLOR,          Hv_gray14,
                NULL);


  FeedbackBSDItem = Hv_VaCreateItem(View,
                Hv_RESIZEPOLICY,       Hv_SHIFTYONRESIZE,
		Hv_TYPE,               Hv_FEEDBACKITEM,
                Hv_BALLOON,            (char *)"Data from mouse position, when at BSD.",
                Hv_PARENT,             Box6,
                Hv_LEFT,               Box6->area->left,
                Hv_TOP,                Box6->area->top,
                Hv_WIDTH,              213,
                Hv_HEIGHT,             100,
                Hv_FEEDBACKDATACOLS,   37,
                Hv_FEEDBACKENTRY,      16, Hv_fixed2, Hv_darkBlue, "",
		Hv_FEEDBACKENTRY,      17, Hv_fixed2, Hv_darkBlue, "",
                Hv_FEEDBACKENTRY,      18, Hv_fixed2, Hv_darkBlue, "",
		Hv_FEEDBACKENTRY,      Hv_SEPARATORTAG,
		Hv_FEEDBACKENTRY,      19, Hv_fixed2, Hv_black, "",
		Hv_FEEDBACKENTRY,      20, Hv_fixed2, Hv_black, "",
                Hv_FEEDBACKENTRY,      21, Hv_fixed2, Hv_black, "",
                Hv_FEEDBACKMARGIN,     20,
                Hv_FILLCOLOR,          Hv_white,
                NULL);


  BSDItem = Hv_VaCreateItem(View,
                Hv_USERDRAW,           BSDUnfolded,
		Hv_BALLOON,            (char *)"The BSD data - not drawn to scale. Solid colors represent ADC hits.  Hatching represents ADC hit within TDC window.",
       		Hv_LEFT,               hotrect->left,
                Hv_TOP,                hotrect->top,
                Hv_WIDTH,              576,
                Hv_HEIGHT,             250,
                NULL);
}

static void   ControlUPVItems(Hv_View View)
{
  Hv_Item	Item;
  Hv_Rect	*hotrect = View->hotrect;
			   

  Button = Hv_VaCreateItem(View,
                Hv_RESIZEPOLICY,        Hv_SHIFTYONRESIZE,
		Hv_TYPE,		Hv_BUTTONITEM,
		Hv_BALLOON,		(char *)"Button to go to the next event.",
		Hv_SINGLECLICK,		ReadForward,
		Hv_LEFT,		hotrect->left + 50,
		Hv_TOP,			hotrect->bottom + 10,
		Hv_WIDTH,		96,
		Hv_HEIGHT,		30,
		Hv_FILLCOLOR,		Hv_white,
		Hv_FONT,		Hv_fixed,
		Hv_TEXT,		"Next event",
		NULL);

  Box6 = Hv_VaCreateItem(View,
                Hv_RESIZEPOLICY,        Hv_SHIFTYONRESIZE,
		Hv_TYPE,		Hv_BOXITEM,
		Hv_LEFT,		hotrect->left,
		Hv_TOP,			hotrect->bottom + 50,
		Hv_WIDTH,		213,
		Hv_HEIGHT,		100,
		Hv_FILLCOLOR,		Hv_gray14,
		NULL);
		

 
  FeedbackUPVItem = Hv_VaCreateItem(View,
		Hv_RESIZEPOLICY,        Hv_SHIFTYONRESIZE,
		Hv_TYPE,		Hv_FEEDBACKITEM,
		Hv_BALLOON,		(char *)"Data from mouse position, when at upv.",
		Hv_PARENT,		Box6,
		Hv_LEFT,		Box6->area->left ,
		Hv_TOP,			Box6->area->top ,
		Hv_WIDTH,               213,
                Hv_HEIGHT,              100, 
                Hv_FEEDBACKDATACOLS,    25,
                Hv_FEEDBACKENTRY,       13, Hv_fixed2, Hv_darkBlue, "", 
                Hv_FEEDBACKENTRY,       14, Hv_fixed2, Hv_black, "",
                Hv_FEEDBACKENTRY,       15, Hv_fixed2,  Hv_black, "", 
                Hv_FEEDBACKENTRY,       16, Hv_fixed2,  Hv_black, "", 
                Hv_FEEDBACKENTRY,       17, Hv_fixed2,  Hv_black, "", 
		Hv_FEEDBACKMARGIN,      20,
                Hv_FILLCOLOR,           Hv_white,
                NULL);

  UPVItem = Hv_VaCreateItem(View,
		Hv_USERDRAW,		UPV,
		Hv_BALLOON,		(char *)"The UPV data- not drawn to scale! Solid colors represent ADC hits. Vertical stripes represent an ADC hit within the TDC window",
		Hv_LEFT,		hotrect->left,
		Hv_TOP,			hotrect->top,
		Hv_WIDTH,		200,
		Hv_HEIGHT,		250,
		NULL);
}

void SetupEdView(Hv_View View)
{
  Hv_Rect	*hotrect = View->hotrect;

  ControlItems(View);
}

void SetupTagView(Hv_View View)
{
  Hv_Rect	*hotrect = View->hotrect;

  ControlTagItems(View);
}

void SetupCPVView(Hv_View View)
{
  Hv_Rect	*hotrect = View->hotrect;

  ControlCPVItems(View);
}

void SetupBSDView(Hv_View View)
{
  Hv_Rect        *hotrect = View->hotrect;

  ControlBSDItems(View);
}

void SetupUPVView(Hv_View View)
{
  Hv_Rect	*hotrect = View->hotrect;

  ControlUPVItems(View);
}
