/*
 * Inputs events from itape and builds a root TTree object
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <disData.h>
#include <dataIO.h>
}

#include <iostream.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <rootclasses.h>

#define ITAPE_BUFSIZE 100000
#define BRANCH_BUFSIZE 100000

VoidFuncPtr_t initfuncs[]={0};

// Initialize the ROOT system
TROOT root("Rint","The ROOT Interactive Interface", initfuncs);

int main(int argc, char *argv[])
{
  char **ifilename;
  char *ofilename;
  int nifiles;

  if (argc<3) {
    printf("usage: %s <input_file_1> ... <output_file>\n",argv[0]);
    printf(" where <input_file_1>... is a list of itape files to be\n");
    printf(" concatenated and saved to root file <output_file>.\n");
    exit(1);
  }
  else {
    ifilename=&argv[1];
    ofilename=argv[argc-1];
    nifiles=argc-2;
  }

  itape_header_t *event=(itape_header_t *)malloc(ITAPE_BUFSIZE);
  if (event == NULL) {
    cerr << "ERROR - failed to allocating input event buffer!" << endl;
    return 1;
  }

  TMCevent *mc_event   =new TMCevent(event,ITAPE_BUFSIZE);
  TLGDhits *lgd_hits   =new TLGDhits(event,ITAPE_BUFSIZE);
  TBSDhits *bsd_hits   =new TBSDhits(event,ITAPE_BUFSIZE);
  TBGVhits *bgv_hits   =new TBGVhits(event,ITAPE_BUFSIZE);
  TCPVhits *cpv_hits   =new TCPVhits(event,ITAPE_BUFSIZE);
  TUPVhits *upv_hits   =new TUPVhits(event,ITAPE_BUFSIZE);
  TTimeList *time_list =new TTimeList(event,ITAPE_BUFSIZE);

  TFile outfile(argv[argc-1],"RECREATE","Radphi data file");
  TTree tree("event_tree","Radphi event tree");
  tree.SetAutoSave();

  TBranch *bMCevent=
     tree.Branch("mc_event","TMCevent",&mc_event,BRANCH_BUFSIZE,0);
  TBranch *bLGDhits=
     tree.Branch("lgd_hits","TLGDhits",&lgd_hits,BRANCH_BUFSIZE,0);
  TBranch *bBSDhits=
     tree.Branch("bsd_hits","TBSDhits",&bsd_hits,BRANCH_BUFSIZE,0);
  TBranch *bBGVhits=
     tree.Branch("bgv_hits","TBGVhits",&bgv_hits,BRANCH_BUFSIZE,0);
  TBranch *bCPVhits=
     tree.Branch("cpv_hits","TCPVhits",&cpv_hits,BRANCH_BUFSIZE,0);
  TBranch *bUPVhits=
     tree.Branch("upv_hits","TUPVhits",&upv_hits,BRANCH_BUFSIZE,0);
  TBranch *bTimeList=
     tree.Branch("time_list","TTimeList",&time_list,BRANCH_BUFSIZE,0);

  for (int ifile=0;ifile<nifiles;ifile++) {
    FILE *fp=fopen(ifilename[ifile],"r");
    if (fp == NULL) {
      cerr << "WARNING - failed to open input events file "
           << ifilename[ifile] << endl;
      tree.Write();
      outfile.Close();
      return 1;
    }
    int i=0;
    while (data_read(fileno(fp),event,ITAPE_BUFSIZE) == DATAIO_OK) {
      cout << ++i << '\r' << flush;
      mc_event->GetGroup();
      lgd_hits->GetGroup();
      bgv_hits->GetGroup();
      bsd_hits->GetGroup();
      cpv_hits->GetGroup();
      upv_hits->GetGroup();
      time_list->GetGroup();
      tree.Fill();
    }
    fclose(fp);
    cout << endl << "copied input file " << ifilename[ifile] << endl;
  }
  outfile.Write();
  outfile.Close();
  return 0;
}
