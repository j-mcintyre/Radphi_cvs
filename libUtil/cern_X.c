#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cernlib.h>

/*
 $Log: cern_X.c,v $
 Revision 1.1  1996/11/21 12:39:18  lfcrob
 Initial revision

*/

static const char rcsid[]="$Id: cern_X.c,v 1.1 1996/11/21 12:39:18 lfcrob Exp $";

#define MAXHISTO 100000

void hplzon_(int *,int *,int *,char *,int);
void hplcon_(int *, int *,int *);
void hplego_(int *,float *,float *);
void iuwk_(int *,int *);
void hplot_(int *,char *,char *,int *,int,int);
void hplint_(int *);

void hplzon(int x, int y, int num, char *opt)
{
  hplzon_(&x, &y, &num, opt, strlen(opt));
}

void hplcon(int histonum, int x, int y)
{
  hplcon_(&histonum, &x, &y);
}

void hplego(int no, float theta, float phi)
{
        hplego_(&no, &theta, &phi);
        return;
}

void iuwk(int num1, int num2)
{
        iuwk_(&num1, &num2);
}

void hplot(int no, char*chopt, char*chcase, int num)
{
        hplot_(&no,chopt,chcase,&num,strlen(chopt),strlen(chcase));
}

void hplint(int no)
{
        hplint_(&no);
}


