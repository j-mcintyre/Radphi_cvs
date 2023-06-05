////////////////////////////////////////////////////////////////////////
//
// main program for BNLGAMS four vector
// dump to file (i.e. command line invocation)
//
// PCF 20 June 1995

////////////////////////////////////////////////////////////////////////

#include "Chees.h"
#include <stdio.h>
#include <stdlib.h>



/* some global variables */
int nAccepted = 0;                 /* counter for accepted events */
int nTotal = 0;                    /* counter for total events */

cheesInfo_t fourVectors;           /* four vecs are stored here per event */
mc_param_t parameters;        /* general input parameters */
particles_db_t *particleDBPtr;         /* particle specifications */



/* some function protypes */
void doStuff(void);
void printInfo(void);
void printHelp(void);


////////////////////////////////////////////////////////////////////////

main(int argc, char *argv[])
{

  char *inputFilename = NULL;
  char *flag = NULL;


// First some command line BS
// Not very interesting, skip down a bit.

  for(int iarg=1; iarg < argc; iarg++)
    if(*argv[iarg] == '-')
      {
        flag = argv[iarg]+1;
	switch(*flag)
          {
	  case 'i':
	    inputFilename = argv[++iarg];
	    break;
	  case 'h':
	    fprintf(stderr, "\n");
	    printHelp();
	    exit(0);
	  default:
	    fprintf(stderr, "\nInvalid command line specification\n");
	    printHelp();
	    exit(EXIT_FAILURE);
	  }
      }

  if(!(inputFilename))  /* bastards didn't give us an input file */ 
    {
      printf("\n\nPlease specify an input file, e.g.,\n\n");
      printHelp();
      exit(EXIT_FAILURE);
    }


  /* here is the crux of the program  - three lines, beautiful */

  CutChees(inputFilename, &parameters, &particleDBPtr);
  doStuff();
  ScrubChees();
  
  /* --------------------------------------------------------- */

  printInfo();

  return 0;
}


////////////////////////////////////////////////////////////////////////


void doStuff(void)
{

  for(int i = 0; i < parameters.codes.nIterations; i++)
    {
      // print out where we are (# of events)
      if(!(i % 50 ))
        {
	  printf("%d\r",i); 
          fflush(stdout);
        }
      
      if(!GenerateOne(&fourVectors))
	++nAccepted;

      ++nTotal;
    }

}


///////////////////////////////////////////////////////////////////////

void printInfo(void)
{

  printf("\n");
  printf("Chees exited successfully\n");
  printf("%d accepted events out of %d total events\n\n", nAccepted, nTotal);

}

////////////////////////////////////////////////////////////////////////

void printHelp(void)
{
  fprintf(stderr, "Usage:  Chees [options]\n");
  fprintf(stderr, "Options are:\n");
  fprintf(stderr, "  -i [filename]\tMakes filename the input param file\n");
  fprintf(stderr, "  -h\t\tPrint this message\n\n");
  fprintf(stderr, "Note: Arguments such as filename, etc. should be ");
  fprintf(stderr, "separateded from the command flags\n\n");
}
