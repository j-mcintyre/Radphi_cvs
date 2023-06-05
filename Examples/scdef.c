#include <string.h>
#include <stdio.h>

#include <cernlib.h>

#include "scalers.h"

void declare_scaler_ntuple (int id, scatable *names, int *p_head, int *p_scalers)
{
  int nvar;
  char title[80];
  char form[12][80];
  float x0[1],x1[1];
  char tags[32];
  int i;

  if (hexist_(&id)) {
    nvar=0;
    hgiven(id,title,&nvar,(char **)form,x0,x1);
    if (strcmp(title,SC_TITLE)) {
      fprintf(stderr,"Warning: possible ntuple version mismatch!\n");
      fprintf(stderr,"title expected: %s\n",SC_TITLE);
      fprintf(stderr,"title in file: %s\n",title);
    }

    hbname(id," ",NULL,"$CLEAR");
    if (p_head != NULL) {
      hbname(id,"head",p_head,"$SET");
    }
    else {
      fprintf(stderr,"declare_ntuple error: null p_head pointer!\n");
      return;
    }
      
    if (p_scalers != NULL) {
      int column=0;
      for (i=0;i<names->nscalers;) {
        int *p_block=&p_scalers[i];
        char columnName[20];
        int last=names->nscalers;
        i=(last-i < 32)? last:i+32;
        sprintf(columnName,"column%d",++column);
        hbname(id,columnName,p_block,"$SET");
      }
    }
    else {
      fprintf(stderr,"declare_ntuple error: null p_scalers pointer!\n");
      return;
    }
  }
  else {
    hbnt(id,SC_TITLE," ");
    if (p_head != NULL) {
      hbname(id,"head",p_head,FORM_HEAD);
    }
    else {
      fprintf(stderr,"declare_ntuple error: null p_head pointer!\n");
      return;
    }
    if (p_scalers != NULL) {
      int column=0;
      for (i=0;i<names->nscalers;) {
        int *p_block=&p_scalers[i];
        char columnName[20];
        char tags[2560];
        int last=names->nscalers;
        last=(last-i < 32)? last:i+32;
        sprintf(tags,"%s:I",names->scaler[i]->tag);
        for (i++;i<last;i++) {
          char tag[80];
          sprintf(tag,",%s:I",names->scaler[i]->tag);
          strcat(tags,tag);
        }
        sprintf(columnName,"column%d",++column);
        hbname(id,columnName,p_block,tags);
      }
    }
    else {
      fprintf(stderr,"declare_ntuple error: null p_scalers pointer!\n");
      return;
    }
  }
}
