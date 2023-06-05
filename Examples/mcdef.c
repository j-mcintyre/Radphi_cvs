#include <stdio.h>
#include <string.h>

#include <mctuple.h>
#include <cernlib.h>

void declare_mc_ntuple(int id, nt_head_t *p_head)
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
      fprintf(stderr,"ERROR: declare_ntuple - p_head pointer was null!\n");
      return;
    }
  }

  else {
    hbnt(id,NT_TITLE," ");
    if (p_head != NULL) {
      hbname(id,"head",p_head,NT_FORM_HEAD);
    }
    else {
      fprintf(stderr,"ERROR: declare_ntuple - p_head pointer was null!\n");
      return;
    }
  }
}

void appendnt_mc_event(int id, nt_mc_event_t *p_mc)
{
  if (hexist_(&id)) {
    if (p_mc != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"mc_event",p_mc,NT_FORM_MC);
      }
      else {
        hbname(id,"mc_event",p_mc,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_mc_event:");
      fprintf(stderr," p_mc pointer was null!\n");
    }
  } 
  else {
    fprintf(stderr,"ERROR: appendnt_mc_event - ntuple id %d does not exist!\n",
            id);
  }
}
