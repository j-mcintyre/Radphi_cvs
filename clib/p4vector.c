/*
  $Log: p4vector.c,v $
  Revision 1.2  1997/05/29 20:26:45  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: p4vector.c,v 1.2 1997/05/29 20:26:45 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void p4vector(const vector4_t *v,FILE* fp)
{
    fprintf(fp,"t: %lf ",v->t);
    p3vector(&v->space,fp);
}
