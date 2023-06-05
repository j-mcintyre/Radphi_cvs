#ifndef TTABLE
#define TTABLE

#include <TROOT.h>

class TString;

class TTable
{
private:
  Int_t fNumColumns;
  Int_t *fColumnWidth;
  Int_t fColumnNum;
  
public:
  TTable(Int_t numColumns, Int_t *columns);
  ~TTable();
  void Print(const char *string);
  void Print(Double_t real);
  void Print(Int_t i);
  void Print(const TString *string);
  void PrintBar();
};

#endif
