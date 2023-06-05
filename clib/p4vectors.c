/*
  $Log: p4vectors.c,v $
  Revision 1.2  1997/05/29 20:26:45  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: p4vectors.c,v 1.2 1997/05/29 20:26:45 radphi Exp $";

#include <stdio.h>
#include <math.h>

#include <clib.h>

void p4vectors(FILE* fp,int npart,const vector4_t *vectors[],const char* title)
{
#define V vectors[i]
	int i;

	float mass;

	fprintf(fp,"\n%s\n\n",title);
	fprintf(fp,"\t N  Energy    Px        Py        Pz        Mass\n\n");
	for (i = 0; i < npart; ++i) {
		mass = v4mag(V);
		fprintf(fp,"\t%2d %9.5f %9.5f %9.5f %9.5f %9.5f\n",i,V->t,V->space.x,V->space.y,
		    V->space.z,mass);
	}
}
