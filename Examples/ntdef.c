#include <string.h>
#include <stdio.h>

#include <ntuple.h>
#include <cernlib.h>

void declare_raw_ntuple(int id, nt_head_t *p_head)
{
  int nvar;
  char title[80];
  char form[12][80];
  float x0[1],x1[1];

  if (hexist_(&id)) {
    nvar=0;
    hgiven(id,title,&nvar,(char **)form,x0,x1);
    if (strcmp(title,NT_TITLE)) {
      fprintf(stderr,"Warning: possible ntuple version mismatch!\n");
      fprintf(stderr,"title expected: %s\n",NT_TITLE);
      fprintf(stderr,"title in file: %s\n",title);
    }
    hbname(id," ",NULL,"$CLEAR");
    if (p_head != NULL) {
      hbname(id,"head",p_head,"$SET");
    }
    else {
      fprintf(stderr,"declare_ntuple error: pointer p_head is null!\n");
      return;
    }
  }

  else {
    hbnt(id,NT_TITLE," ");
    if (p_head != NULL) {
      hbname(id,"head",p_head,NT_FORM_HEAD);
    }
    else {
      fprintf(stderr,"declare_ntuple error: pointer p_head is null!\n");
      return;
    }
  }
}

void appendnt_raw_event(int id, nt_raw_event_t *p_raw)
{

  if (hexist_(&id)) {
    if (p_raw != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"raw",p_raw,NT_FORM_RAW);
      }
      else {
        hbname(id,"raw",p_raw,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_raw_event:");
      fprintf(stderr," p_raw pointer was null!\n");
    }
  } 
  else {
    fprintf(stderr,"ERROR: appendnt_raw_event - ntuple id %d does not exist!\n",
              id);
  }
}
