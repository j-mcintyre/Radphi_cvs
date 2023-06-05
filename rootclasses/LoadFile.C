TFile *file;
TTree *tree;
TMCdata *mcData;
THitsData *hitsData;
TReconstruction *phots;

Bool_t LoadFile(const char *filename)
{
  char treeName[]="h10";
  file=new TFile(filename,"READ",filename);
  tree=(TTree *)file->Get(treeName);
  if (tree!=NULL) {
    mcData=new TMCdata(tree);
    hitsData=new THitsData(tree);
    phots=new TReconstruction(tree);
    return(kTRUE);
  }
  else {
    cout << "Fatal: could not open tree " << treeName << " in file " <<
            filename << endl;
    return(kFALSE);
  }
}
