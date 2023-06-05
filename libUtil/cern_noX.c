#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cernlib.h>
/*
 $Log: cern_noX.c,v $
 Revision 1.4  2001/04/13 16:28:30  radphi
 Changes made by crsteffe@mantrid00
 added hrget definition _craig

 Revision 1.3  1998/06/22 22:20:15  radphi
 Changes made by jonesrt@jlabs4
 fixed a bug in the c routine hgiven that front-ends the hbook call hgiven - RTJ

 * Revision 1.2  1997/05/03  20:38:08  radphi
 * Changes made by lfcrob@dustbunny
 * Removed pawc declaration.
 * Removed setup_histo(...) stuff, replaced with hlimit()
 * as it should be
 *
 * Revision 1.1  1996/11/21  12:39:18  lfcrob
 * Initial revision
 *
*/

static const char rcsid[]="$Id: cern_noX.c,v 1.4 2001/04/13 16:28:30 radphi Exp $";


void hbook1_(int *,char *,int *,float *,float*,float *,int);
void hbook2_(int *,char *,int *,float *,float*,int *,float *,float *,float *,int);
void hbookn_(int *,char *,int *,char *,int *,char *,int,int,int);
void hf1_(int *,float *,float *);
void hf2_(int *,float *,float *,float *);
void hfill_(int *,float *,float *,float *);
void hfn_(int *,float *);
void hidopt_(int *,char *,int);
void hrput_(int *,char *,char *,int,int);
void hrget_(int *,char *,char *,int,int);
void hunpak_(int *,float *,char *,int *,int);
void hpak_(int *,float *);
void hldir_(char *,char *,int,int);
void hmdir_(char *,char *,int,int);
void hcdir_(char *,char *,int,int);
void hropen_(int *,char *,char *,char *,int *,int *,int,int,int);
void hrout_(int *,int *,char *,int);
void hrin_(int *, int *, int *);
void hrend_(char *,int);
void hreset_(int *,char *,int);
void hlimit_(int *);
void hlimap_(int *,char *,int);
void hbnt_(int *id ,char* name, char* chmod,int ,int);
void hbname_(int *id, char* chblok, void*variable, char*chform, int,int);
void hfnt_(int*id);
void hfntb_(int*id,char* chblok, int);
void hprnt_(int*id);
void hgiven_(int *,char *,int *,char *,float *,float *,int,int);
void hgn_(int *,int *,int *,float *,int *);
void hgnpar_(int *,char *,int);
void hgnf_(int *,int *,float *,int *);
void hnoent_(int *,int *);
void hgnt_(int *,int *,int *);
void hgntb_(int *,char *,int *,int *,int);
void hgntv_(int *,char *,int *,int *,int *,int);
void hgntf_(int *,int *,int *);
void hprntu_(int *);


void hlimit(int size)
{
  hlimit_(&size);
  return;
}
/*
void hlimap(int size, char *name)
{
  hlimap_(&size,name,strlen(name));
  return;
}
*/

void hbookn(int id,char *title,int ntags,char *dir,int block_size,char *tags[])
{ 
  int i;
  int j;
  char *ftag = NULL;
  ftag = malloc(ntags *8* sizeof(char) +1 );
  for(i=0;i<ntags;i++)
    {
      memcpy(ftag+i*8,tags[i],strlen(tags[i]));
    }
  hbookn_(&id,title,&ntags,dir,&block_size,ftag,strlen(title),strlen(dir),8L);
}
  
  

void hf1(int id,float data,float weight)
{
  hf1_(&id,&data,&weight);
}


void hf2(int id,float data1,float data2,float weight)
{
  hf2_(&id,&data1,&data2,&weight);
}

void hfn(int id, float data[])
{
  hfn_(&id,&data[0]); 
}

void hidopt(int id,char *command)
{
  hidopt_(&id,command,strlen(command));
  return;
}

void hrput(int id, char *fname,char *command)
{
  hrput_(&id,fname,command,strlen(fname),strlen(command));
}

/* hrget added 2001 April 13 by Craig Steffen */
void hrget(int id, char *fname,char *command)
{
  hrget_(&id,fname,command,strlen(fname),strlen(command));
}

void save_histo(int id, char fname[80])
{
  /* R. Lindenbusch's saver */
  char command[80];
  char command2[80];
  
  sprintf(command,"BLAC");
  hidopt_(&id,command,strlen(command));
  sprintf(command,"1EVL");
  hidopt_(&id,command,strlen(command));
  sprintf(command,"NT");
  hrput_(&id,fname,command,strlen(fname),strlen(command));
}

void hunpak(int histo,float contents[],char choice[], int num)
{
  hunpak_(&histo,contents,choice,&num,strlen(choice));
}

void hpak(int histo,float contents[])
{
  hpak_(&histo,contents);
}

void hldir(char dir[],char flag[])
{
  hldir_(dir,flag,strlen(dir),strlen(flag));
}

void hmdir(char dir[],char flag[])
{
  hmdir_(dir,flag,strlen(dir),strlen(flag));
}

void hcdir(char dir[],char flag[])
{
  hcdir_(dir,flag,strlen(dir),strlen(flag));
}

void hropen(int lun, char *name, char*filename, char*status, int *lrec, int *istat)
{
  hropen_(&lun, name, filename, status, lrec, istat, strlen(name),
	  strlen(filename), strlen(status));
}

void hrout(int num, int icycle, char*opt)
{
  hrout_(&num, &icycle, opt, strlen(opt));
  return;
}

void hrin(int id, int icycle, int iofset)
{
  hrin_(&id, &icycle, &iofset);
}

void hrend(char*filename)
{
  hrend_(filename, strlen(filename));
  return;
}

void hreset(int no, char* opt)
{
  hreset_(&no, opt, strlen(opt));
  return;
}

void hbook2(int no, char*name, int xbins, float xmin, float xmax, int ybins,
	    float ymin, float ymax, float weight)
{
  hbook2_(&no, name, &xbins, &xmin, &xmax, &ybins,&ymin, &ymax, &weight, strlen(name));
}


void hbook1(int no, char*name, int nbins, float min, float max, float v)
{
  hbook1_(&no, name, &nbins, &min, &max, &v, strlen(name));
}

void hfill(int no, float xvalue, float yvalue, float weight)
{
  hfill_(&no, &xvalue, &yvalue, &weight);
}

void hbnt(int id,char*chtitle,char*chopt)
{
  hbnt_(&id,chtitle,chopt,strlen(chtitle),strlen(chopt));
}

void hbname(int id,char*chblok,void*variable,char*chform)
{
  hbname_(&id,chblok,variable,chform,strlen(chblok),strlen(chform)); 
}

void hfnt(int id)
{
  hfnt_(&id);
}
 
void hfntb(int id,char* chblok)
{
  hfntb_(&id,chblok,strlen(chblok));
}

void hprnt(int id)
{
  hprnt_(&id);
}

/* Note: because nvar is to be returned you need to pass the address into
   hgiven! */
void hgiven(int id, char *chtitl, int *nvar, char *chtag[], float rlow[], float rhigh[])
{
  int i;
  char *ftag = NULL;
  char title[1024];
  int num;

  num = *nvar;
  ftag = malloc(num *8* sizeof(char) +1 );
  hgiven_(&id,title,&num,ftag,&rlow[0],&rhigh[0],1024,8L);
  for(i=0;i<num;i++)
    if(i<*nvar)
      memcpy(chtag[i],ftag+i*8,8);
  *nvar = num;
  /* Since fortran fills the rest of the title with spaces and will not let me
     find out how long the title string is, I need to guess at the end of the
     title and cut the rest off --> FORTRAN SUCKS!!! */
  i=0;
  while(strncmp(&title[i],"        ",8))
    i++;
  strncpy(chtitl,title,i);
  chtitl[i]=0;
}

void hgn(int id, int idn, int idnevt, float x[], int *ierror)
{
  hgn_(&id,&idn,&idnevt,&x[0],ierror);
}

void hgnpar(int id, char *chrout)
{
  hgnpar_(&id,chrout,strlen(chrout));
}

void hgnf(int id, int idnevt, float x[], int *ierror)
{
  hgnf_(&id,&idnevt,&x[0],ierror);
}

void hnoent(int id, int *noent)
{
  hnoent_(&id,noent);
}

void hgnt(int id, int idnevt, int *ierr)
{
  hgnt_(&id,&idnevt,ierr);
}

void hgntb(int id, char *chblok, int irow, int *ierr)
{
  hgntb_(&id,chblok,&irow,ierr,strlen(chblok));
}

void hgntv(int id, char *chvar[], int nvar, int irow, int *ierr)
{
  int i;
  int offset = 0; 
  char *ftag = NULL;

  ftag = malloc(sizeof(chvar)+1);
  for(i=0;1<nvar;i++) {
    memcpy(ftag+offset,chvar[i],strlen(chvar[i]));
    offset += (int) strlen(chvar[i]);
  }
  hgntv_(&id,ftag,&nvar,&irow,ierr,strlen(ftag));
}

void hgntf(int id, int irow, int *ierr)
{
  hgntf_(&id,&irow,ierr);
}

/* Function to print info about ntuple. */
void hprntu(int id)
{
  hprntu_(&id);
}
