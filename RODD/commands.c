/* commands.c */
/*
 * This is the bulk of the command line module. Commands and their
 * actions are defined here. 
 * Command Line is the loop we continually execute until told to go
 * or to quit.
 *
 */

/*  DSA - added EXIT and E as synonyms for QUIT     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <time.h>

#define LINE_SIZE 256	    /* size of the line we get from user */
#define NUM_FILES 5	    /* number of files deep we are allowed to go */

#define UNKNOWN_COMMAND -1
#define SET 0
#define LOAD 1
#define ADDRUN 2
#define DELRUN 3
#define LIST 4
#define CLEAR 5
#define GO 6
#define QUIT 7
#define Q 8
#define HELP 9
#define EXIT 10
#define E 11
#define VERBOSE 12
#define QUIET 13
#define NUM_COMMANDS 14

#define LOG 0
#define INPUT 1
#define OUTPUT 2
#define TRIGGER 3
#define HIST 4
#define NUM_OBJECTS 5
#define NUM_TRIGGERS 6
#define UNKNOWN_OBJECT -1

#define PHYSICS 0
#define LASER 1
#define BGV  2
#define PEDESTAL 3
#define BASETEST 4
#define ALL 5


#include "runs.h"
#include "errors.h"
#include "global.h"

FILE* infile; /* = stdin;	*/ /* dsa test */           /* current input file pointer */
static FILE *file_stack[NUM_FILES]; /* stack of file pointers */
static int top = -1;		    /* top of our stack */

static char *commands[NUM_COMMANDS] = /* strings of commands known to command line */
        {"set", "load", "addrun", "delrun", 
        "list", "clear", "go", "quit", "q", "help", "exit", "e", "verbose", "quiet"};

static char *set_objects[NUM_OBJECTS] = {"log", "input", "output", "trigger", "hist"};  /* objects that can be "set" */

static char *set_triggers[NUM_TRIGGERS] = {"physics", "laser", "bgv", "pedestal", "basetest", "all"};  /* triggersthat can be "selected" */

static char general_help[] = {"General Help message\nTo get help with\
specific commands type: help <command>\n"};	

static char *specific_help[NUM_COMMANDS] = {
   	"SET---\nSyntax: set <object> <value>\n\
Synopsis: Changes the value of an internal object\n\
          Legal objects: log, input, output, hist, trigger\n", 
	
	"LOAD\n----\nSyntax: load <filename>\n\
Synopsis: Reads in and executes commands located in <filename>", 
	
	"ADDRUN\n------\nSyntax: addrun <runnumber>\
{<runnumber>}\n\taddrun <runnumber>\
(<eventnum>-<eventnum>)\n\taddrun <runnumber>\
(<eventnum>)\nSynopsis: Puts the specified\
 run number and\n\toptional event cuts on the list to be processed", 
    	
	"DELRUN\n------\nSyntax: delrun <runnmuber>[<runnumber>]\n\
Synopsis: Removes all entries of the specified run numbers from list", 
	
	"LIST\n----\nLists all the jobs to be processed", 
    	
	"CLEAR\n-----\nRemoves all jobs on list", 
	
	"GO\n--\nBegin processing runs on list.", 
	
	"QUIT\n----\nExit to OS", 
	
	"Q\n-\nSame as QUIT", 
	
	"HELP\n----\nLists all legal commands",     	

        "EXIT\n----\nSame as QUIT",     

        "E\n----\nSame as QUIT",     

        "VERBOSE\n----\nWrite unpacking info to LOG file WARNING - this gets lengthy very quickly",     

        "QUIET\n----\nTurns off VERBOSE mode",     

};

static void
list_commands(){
  int i;
  fprintf(stdout, "The current legal commands understood by the system are:\n");
  for (i = 0; i < NUM_COMMANDS; i++) {
    fprintf(stdout, "%s ", commands[i]);
  }
}

static void 
generic_help () {
  fprintf(stdout, "\n%s",  general_help);
  list_commands();	
}

static void 
help_on( int the_command ) {
  if (the_command < 0) {
    fprintf(stdout, "Cannot get help on that unknown command.\n");
    list_commands();
  }
  else fprintf(stdout, "\n%s",  specific_help[the_command]);
}

static int 
process_command(char theinput[]) { 
  int i;
  for (i = 0; i < NUM_COMMANDS ; i++)
    if (!strcmp(theinput, commands[i])) 
      return(i);    
  return(UNKNOWN_COMMAND);
}

static void 
addruns(char list[]) {
    int in_paren = 0;
    int dash = 0;
    int number_read = 0;
    int star = 0;
    char curr[11];
    int runnum;
    long lowbnd;
    char c;
    int j = 0;
    int k = 0;
    
    strcpy(curr, "\0");
    c = list[k++];
    while (c != '\0') {
	if (isdigit(c)) {	    /* a number */
	    curr[j] = c;
	    j += 1;
	    number_read = 1;
	}
	else {			    /* a character */
	    switch (c) {
		case '(' :	    /* open paren */
		    if (!number_read) goto adderror; /* no run given */
		    curr[j] = '\0';
		    runnum = atoi(curr);   /* run read, now event cuts */
		    j = 0;
		    in_paren = 1; number_read = 0; dash = 0; star = 0;
		    break;  
		
		case '-' :		/* something thru something */
		    if (!number_read) goto adderror; /* no number */
		    if (!in_paren) {	/* runnum - runnum */
			curr[j] = '\0';
			runnum = atoi(curr);
		    }
		    else {		/* runnum(event-event) */
			curr[j] = '\0';
			lowbnd = atoi(curr);
		    }
		    j = 0;
		    dash = 1; number_read = 0;
		    break;
		
		case ')' :	    /* end of run cut */
		    if (star) {
			star = 0; in_paren = 0; break;
		    }
		    if (!in_paren || !number_read) 
			goto adderror;	/* mis-matched parens */
		    if (dash) {	    /* runnum(event - event) */
			curr[j] = '\0';
			if (lowbnd > atoi(curr))
			    size_error(lowbnd, atoi(curr));
			else
			    set_multi_events(runnum, lowbnd, atoi(curr));
		    }
		    else {	    /* runnum(event) */
			curr[j] = '\0';
			set_event(runnum, atoi(curr));
		    }
		    j = 0;
		    dash = 0; number_read = 0; in_paren = 0;
		    break;
		
		case ' ' :	    /* whitespace */
		    if (number_read) {	/* end of a number */
			if (in_paren) { 
			    if (dash) {	/* runnum(event-event  */
				curr[j] = '\0';
				if (lowbnd > atoi(curr))
				    size_error(lowbnd, atoi(curr));
				else
				    set_multi_events(runnum, lowbnd, atoi(curr));
				dash = 0;
			    }
			    else {     /* runnum(event */
				curr[j] = '\0';
				set_event(runnum, atoi(curr));
			    }
			}
			else {		/* runnum-runnum */
			    if (dash) {
				curr[j] = '\0';
				set_multi_runs(runnum, atoi(curr));
				dash = 0;
			    }
			    else {     /* runnum */
				curr[j] = '\0';
				runnum = atoi(curr);
				set_run(runnum);
			    }
			}
		    }
		    star = 0; number_read = 0; j = 0;
		    break;
		
		case '*' :	    /* go till end of run */
		    if (!in_paren && !dash) goto adderror;
		    else {
			set_multi_events(runnum, lowbnd, INT_MAX);
		    }
		    star = 1; number_read = 0; j = 0; dash = 0;
		    
		default :  goto adderror;
	    }			/* switch */
	}			/* else */
	c = list[k++];
    }				/* while */
    
    if (in_paren) { /* assume they meant to close the paren but forgot */
	if (dash) {
	    curr[j] = '\0';
	    if (lowbnd > atoi(curr))
		size_error(lowbnd, atoi(curr));
	    else
		set_multi_events(runnum, lowbnd, atoi(curr));
	}
	else {
	    curr[j] = '\0';
	    set_event(runnum, atoi(curr));
	}
    }
    else if (number_read) { /* a number was the last thing read */
	if (dash) {
	    curr[j] = '\0';
	    if (runnum > atoi(curr))
		size_error(runnum, atoi(curr));
	    else
		set_multi_runs(runnum, atoi(curr));
	}    
	else {
	    curr[j] = '\0';
	    runnum = atoi(curr);
	    set_run(runnum);
	}
    }	
    return;
    adderror: {
	addrun_error(c);
	return;
    }		    /* error */
}		    /* addruns */

static void 
delruns(char items[]){
    char number[11];
    int lower;
    int number_read = 0;
    int dash = 0;
    int size;
    char c;
    int i = 0;
    int j = 0;
    c = items[i++];
    strcpy(number, "\0");
    while (c != '\0') {
	if (isdigit(c)) {
	    number[j++] = c;
	    number_read = 1;
	} /* if isdigit */
	else {
	    switch (c) {
		case '-' :
		   if (!number_read) goto delerror;
		   number[j] = '\0';
		   lower = atoi(number);
		   dash = 1; number_read = 0; j = 0;
		   break;
		   
		case ' ' :
		    if(dash) {
			number[j] = '\0';
			if (lower > atoi(number))
			    size_error(lower, atoi(number));
			else
			    multi_remove(lower, atoi(number));
		    }
		    else {
			number[j] = '\0';
			remove_run(atoi(number));
		    }
		    dash = 0; number_read = 0; j = 0;
		    break;
		
		default :
		    delrun_error(c);
		    return;
	    } /* switch */
	} /* else */
	c = items[i++];
    } /* while */
    if (number_read) {
	if (dash) {
	    number[j] = '\0';
	    if (lower > atoi(number))
		size_error(lower, atoi(number));
	    else
		multi_remove(lower, atoi(number));
	}
	else {
	    number[j] = '\0';
	    remove_run(atoi(number));
	}
    } /* if number read */
    return;
    delerror : {
	delrun_error(c);
	return;
    }
} /* delruns */

static void 
import_file( char filename[] ) {
  if ((top + 1) == NUM_FILES) stack_error(); /* error, too many files */
  else {
    /*    infile = stdin; */
    file_stack[++top] = infile;
    if ((infile = fopen(filename, "r")) == NULL)
      file_error(filename);
    BATCH = 1;
    fprintf(logfi, "Importing command file: %s\n", filename);
  }
}

static int
process_object(char *obj) {
  int i;
  for (i = 0; i < NUM_OBJECTS; i++) {
    if (!strcmp(obj, set_objects[i]))
      return(i);
  }
  return(UNKNOWN_OBJECT);
}
static int
process_trigger(char *obj) {
  int i;
  for (i = 0; i < NUM_TRIGGERS; i++) {
    if (!strcmp(obj, set_triggers[i]))
      return(i);
  }
  return(UNKNOWN_OBJECT);
}

static void 
set_command(char args[]) {
  char object[20], value[100];
  FILE* temp;
  time_t curr;
  int obj;
  int trig;

  if (sscanf(args, "%s %s", object, value) < 2) {
    set_error(); 
  }

  obj = process_object(object);

  switch(obj) {
  case LOG:
    fprintf(logfi, "Log file closed and switched to %s\n", value);
    if ((temp = fopen(value, "a")) == NULL)
      file_error(value);
    fclose(logfi);
    logfi = temp; temp = NULL;
    time(&curr);
    fprintf(logfi, "New logfile started at %s", ctime(&curr));
    if (!BATCH) {
      fprintf(stdout, "Logfile switched to %s\n", value);
    }
    break;
  case INPUT:
    if (input_tape)
      free(input_tape);
    input_tape = malloc(strlen(value)+1);
    strcpy(input_tape, value);
    fprintf(logfi, "Input tape name changed to: %s\n", input_tape);
    if (!BATCH)
      fprintf(stdout, "Input tape name set to: %s\n", input_tape);
    break;
  case OUTPUT:
    if (output_tape)
      free(output_tape);
    output_tape = malloc(strlen(value)+1);
    strcpy(output_tape, value);
    fprintf(logfi, "Output tape name changed to: %s\n", output_tape);
    if (!BATCH)
      fprintf(stdout, "Output tape name set to: %s\n", output_tape);
    break;
  case HIST:
    if (histo_file_name)
      free(histo_file_name);
    histo_file_name = malloc(strlen(value)+1);
    strcpy(histo_file_name, value);
    fprintf(logfi, "Histogram name changed to: %s\n", histo_file_name);
    if (!BATCH)
      fprintf(stdout, "Histogram name  name set to: %s\n", histo_file_name);
    break;
  case TRIGGER:
    trig = process_trigger(value);
    switch(trig) {
      case PHYSICS:
       trig_flag=2;
       break;
      case LASER:
       trig_flag=4;
       break;
      case BGV:
       trig_flag=8;
       break;
      case PEDESTAL:
       trig_flag=16;
       break;
      case BASETEST:
       trig_flag=32;
       break;
      case ALL:
       trig_flag=0;
      default:
       trig_flag=0;
    }
    fprintf(logfi, "Only histogram triggers of type: %s\n", value);
    if (!BATCH)
      fprintf(stdout, "Only histogram triggers of type: %s\n", value);
    break;
  default:
    set_error();
    
  }
}

static void 
execute_command(int to_do, char *argument) {
    switch (to_do) {
    case ADDRUN : addruns(argument); break;
    case Q :
    case E :
    case EXIT :
    case QUIT : break;
    case CLEAR : clear_runs(); break;
    case DELRUN : delruns(argument); break;
    case GO : fprintf(logfi, "\n****** Beginning data analysis. ******\n"); show_listing(logfi); break;
    case HELP : 
      if (strlen(argument) == 0) generic_help();
      else {
	int request;
	char keyword[25];
	sscanf(argument, "%s", keyword);
	request = process_command(keyword);
	help_on(request);
      }
      break;
    case LIST : 
      if (BATCH) show_listing(stderr);
      else show_listing(stdout);
      break;
    case VERBOSE:
      fprintf(logfi, "Verbose mode on\n");
      unpackSetVerbose(1);
      break;
    case QUIET:
      fprintf(logfi, "Verbose mode off\n");
      unpackSetVerbose(0);
      break;
    case LOAD : import_file(argument); break;
    case SET : set_command(argument); break;
    case UNKNOWN_COMMAND : unknown_command(); break;
    default : fprintf(stderr, "Error in command parser\n");
      fprintf(stderr, "RODD bailing out.\n");
      fprintf(logfi, "RODD terminating due to error in command parser.\n");
      fclose(logfi);
      exit(EXIT_FAILURE);
    }
}

int 
command_line() {
  char input[LINE_SIZE];
  char word[100];
  int action;
  while(1) {
    if (feof(infile)) { 
      if (infile == stdin) return 0;
      fclose(infile);
      if (top == -1) infile = stdin;
      else infile = file_stack[top--]; 
      if (infile == stdin) BATCH = 0;
    }
    if (infile == stdin) printf("\nRODD%%: ");
    if ((fscanf(infile, "%[^\n]", input)) <= 0) {
      fgetc(infile);
      continue;
    }
    else fgetc(infile);         /* get the \n outta the way */
    sscanf(input, "%s", word);
    action = process_command(word);
    strcpy(input, &input[strlen(word) + 1]);
    execute_command(action, input);
    if ((action == QUIT) || (action == Q) || (action == EXIT) || (action == E)){
      return 0;
    }
    if (action == GO) {
      return 1;
    }
  } /* end while(1) */
}
