//--------------------------------------------------
#include "TROOT.h"
#include "TRint.h"

#include "TTimeList.h"
#include "TBSDhits.h"
#include "TCPVhits.h"
#include "TUPVhits.h"
#include "TBGVhits.h"
#include "TLGDhits.h"

int Error; //left undefined by Motif

extern void  InitGui();  // initializer for GUI needed for interactive interface
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };

// Initialize the ROOT system
TROOT root("Rint","The ROOT Interactive Interface", initfuncs);

int main(int argc, char **argv)
{
    // Create interactive interface
    TRint *theApp = new TRint("ROOT example", &argc, argv, NULL, 0);
    TTimeList t;

    // Run interactive interface
    theApp->Run();

    return(0);
}
//--------------------------------------------------
