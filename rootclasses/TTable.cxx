#include <iostream.h>
#include <stdio.h>
#include <TString.h>
#include "TTable.h"

TTable::TTable(Int_t numColumns,Int_t *columns)
	:fNumColumns(numColumns),fColumnNum(0)
{
  fColumnWidth=new Int_t[numColumns];
  memcpy(fColumnWidth,columns,sizeof(Int_t)*numColumns);
}

TTable::~TTable()
{
  delete [] fColumnWidth;
}

void TTable::Print(const char *string)
{
  TString s=string;
  if (s.Length()>fColumnWidth[fColumnNum]) 
    s=s(0,fColumnWidth[fColumnNum]);
  while(s.Length()<fColumnWidth[fColumnNum]) s+=" ";
  s+=" | ";
  cout << s;
  if (++fColumnNum>=fNumColumns) {
    cout << endl;
    fColumnNum=0;
  }
}

void TTable::Print(Double_t real)
{
  char string[30];
  sprintf(string,"%lf",real);
  Print(string);
}

void TTable::Print(Int_t i)
{
  char string[50];
  sprintf(string,"%d",(int)i);
  Print(string);
}

void TTable::Print(const TString *string)
{
  Print(string->Data());
}

void TTable::PrintBar()
{
  for(Int_t columnNum=0;columnNum<fNumColumns;columnNum++) {
    for(Int_t i=0;i<fColumnWidth[columnNum];i++)
      cout << '-';
    cout << "-+-";
  }
  cout << endl;
}

#ifdef VERIFIER

main()
{
  Int_t cw[]= {4,5,6};
  TTable table(3,cw);
  table.Print(0);
  table.Print(0.0);
  table.Print("hello");
  table.Print("goodbye");
  table.Print(53443);
  table.Print(5.33433);
  table.Print(1.0/3.0);
  table.Print("dd");
  table.Print("");
}

#endif
