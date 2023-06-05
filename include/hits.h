#ifndef HITS_H_INCLUDED
#define HITS_H_INCLUDED 

#include <ntypes.h>
#include <disData.h>

#define BUFSIZE 100000
int makeHits(itape_header_t *event); 
int configureIt(itape_header_t *event); 

#endif
