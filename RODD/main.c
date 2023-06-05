/* main.c */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "global.h"
#include "errors.h"
#include <string.h>
#include <disIO.h>

extern void setupHbook();
extern void closeHbook();

/* %%% all located in global.h %%% */
FILE* logfi = NULL;
int BATCH = 0;
char *input_tape = 0;
char *output_tape = 0;
char *dispatcher = 0;
char *comfile = 0;
char *logfile = 0;
char *histo_file_name = 0;

int trig_flag;


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

extern void analysis();
extern int command_line(void);

static time_t now;

static char default_log[] = "rodd.log";

static void 
startup() {

  setupHbook(); 

  if (logfile != NULL) 
    logfi = fopen(logfile, "a");
  else
    logfi = fopen(default_log, "a");
  if (logfi)
    fprintf(logfi, "RODD started up for %s on %s.\n", getenv("USER"), 
            getenv("HOST"));
  else {
    fprintf(stderr, "Error in startup: unable to open log file.\n");
    exit(9);
  }
  now = time(NULL);
  fprintf(logfi, "Startup time is: %s", ctime(&now));	

  infile = stdin;
  histo_file_name = malloc(999);
  strcpy(histo_file_name,"radphi.hbook");

  if (comfile) {
    if ((infile = fopen(comfile, "r")) == NULL)
      file_error(comfile);
  }
  
  if (dispatcher) fprintf(logfi,"Connecting to Dispatcher %s\n",dispatcher);
  if (input_tape) fprintf(logfi, "ITAPE input: %s\n", input_tape);
  if (output_tape) fprintf(logfi, "ITAPE output: %s\n", output_tape);
  if (comfile) fprintf(logfi, "Commands file: %s\n", comfile);
}

static void 
shutdown() {
  now = time(NULL);
  fprintf(logfi, "RODD finished at: %s", ctime(&now));
  fprintf(logfi, "----------------------------------------------\n");
  fclose(logfi);
  closeHbook();  
}

static void
Usage(char* argv0) {
  fprintf(stderr, "%s <options>\n", argv0);
  fprintf(stderr, "Options: \n");
  fprintf(stderr,"\t-d Dispatcher\n");
  fprintf(stderr, "\t-i itapeInputFile\n");
  fprintf(stderr, "\t-o itapeOutputFile\n");
  fprintf(stderr, "\t-c batchCommandFile\n");
  fprintf(stderr, "\t-l logFilename\n");
  exit(1);
}

int 
main(int argc, char** argv) {
  
  char *argptr;
  int i;
  for (i = 1; i < argc; i++) {
    if (*argv[i] == '-') {
      argptr = argv[i] + 1;
      switch(*argptr) {
      case 'd':
	if (argv[++i] == NULL)
           Usage(argv[0]);
        if (dispatcher)
           free(dispatcher);
        dispatcher = malloc(strlen(argv[i])+1);
        strcpy(dispatcher,argv[i]);
	break;	
      case 'i':
	if (argv[++i] == NULL)
           Usage(argv[0]);
        if (input_tape)
           free(input_tape);
        input_tape = malloc(strlen(argv[i])+1);
        strcpy(input_tape,argv[i]);
	break;
      case 'o':
	if (argv[++i] == NULL)
           Usage(argv[0]);
        if (output_tape)
           free(output_tape);
        output_tape = malloc(strlen(argv[i])+1);
        strcpy(output_tape,argv[i]);
	break;
      case 'c':
	if (argv[++i] == NULL)
           Usage(argv[0]);
        if (comfile)
           free(comfile);
        comfile = malloc(strlen(argv[i])+1);
        strcpy(comfile,argv[i]);
	BATCH = 1;
	break;
      case 'l':
	if (argv[++i] == NULL)
           Usage(argv[0]);
        if (logfile)
           free(logfile);
        logfile = malloc(strlen(argv[i])+1);
        strcpy(logfile,argv[i]);
	break;
      case 'h':
	Usage(argv[0]);
	break;
      default:
	fprintf(stderr, "Unknown argument : %s\n", argptr);
	Usage(argv[0]);
      } /* end switch */
    }
    else {
      fprintf(stderr, "Unknown argument : %s\n", argv[i]);
      Usage(argv[0]);
    }
  }
  if(dispatcher && (input_tape || output_tape)){
    fprintf(stderr,"Input/Ouput tapes not allowed in dispatcher mode.\n");
    Usage(argv[0]);
  }
  if(!dispatcher){
    if (((input_tape == NULL) || (output_tape == NULL)) && (comfile == NULL)) {
      fprintf(stderr, "If ITAPE going to stdout, or coming from stdin, must specify command file\n");
      Usage(argv[0]);
    }
  }

  /* begin real work */

  startup();
  
  /* the main driving loop */
  while(command_line()) 
    analysis();
  shutdown();
  /* cleanup */
  return(EXIT_SUCCESS);
}
