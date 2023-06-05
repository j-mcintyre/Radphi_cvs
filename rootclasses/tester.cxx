//--------------------------------------------------
#include <TROOT.h>
#include <TRint.h>

#include "TMCdata.h"
#include "THitsData.h"
#include "TReconstruction.h"
#include "TChannel.h"
#include "TResolver.h"

int main(int argc, char **argv)
{
    // Initialize the ROOT system
    TROOT root("simple","ROOT in batch mode");

    TFile file("work/simphi3eta.root","READ","MC test data");
    TTree *tree=(TTree *)file.Get("h10");
    TReconstruction phots(tree);
    phots.FindSigma();
    
    file.Close();

    return(0);
}
//--------------------------------------------------
