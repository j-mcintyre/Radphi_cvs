/* errors.c */

#include <stdio.h>
#include <stdlib.h>
#include "errors.h"
#include "global.h"
#include "runs.h"

void umapload_error(int run) {
  fprintf(stderr, "Umap failed to load db for run: %d\n", run);
  fprintf(stderr, "Aborting remainder of decoding\n");
  fprintf(logfi, "ERROR: Umap failed to load db for run :%d\n", run);
  fprintf(logfi, ">>>>>> Not attemtpting to continue with processing\n");
}

void umap_error(int ROC, int slot, int channel) {
  fprintf(stderr, "Error in getting hardware mapping information, check log\n");
  fprintf(logfi, "ERROR: The following did not map correctly: \n");
  fprintf(logfi, "   ROC: %d, Slot: %d, Channel %d\n", ROC, slot, channel);
  fprintf(logfi, "   Either incorrect slot to module mapping, or entry does not exist\n");
}

void malloc_error( char* source, int line ) {
    fprintf(stderr, "\n\nError in memory allocation!\n");
    fprintf(stderr, "Memory allocation failed at line %d in source file %s\n", line, source);
    fprintf(stderr, "Unable to continue, bailing out.\n");
    fprintf(logfi, "RODD terminated due to malloc error\n");
    fprintf(logfi, "----------------------------------------------\n");
    fclose(logfi);
    exit(EXIT_FAILURE);
}

void stack_error() {
    fprintf(stderr, "\n\nError in depth of command files.\n");
    fprintf(stderr, "Recursive depth of command files too deep.\n");
    fprintf(stderr, "Edit commands.c to change depth.\n");
    fprintf(stderr, "Terminating RODD.\n");
    fprintf(logfi, "RODD terminated due to stack error\n");
    fprintf(logfi, "----------------------------------------------\n");
    fclose(logfi);
    exit(EXIT_FAILURE);
}

void delrun_error( char c ) {
    fprintf(logfi, "Error found in delrun command with character %c\n", c);
    switch (c) {
	case '-': fprintf(stderr, "Missing a number before dash!\n"); break;
	default : fprintf(stderr, "Illegal character %c.\n", c);
    }
    if (!BATCH) {
	fprintf(stderr, "Rest of line ignored\n");
	return;
    }
    else {
	fprintf(stderr, "\n\nError in delrun command in command file!\n");
	fprintf(stderr, "Fix syntax of command file\n");
	fprintf(stderr, "RODD terminating.");
	fprintf(logfi, "RODD terminated due to delrun syntax error in file.\n");
	fprintf(logfi, "----------------------------------------------\n");
	fclose(logfi);
	exit(EXIT_FAILURE);
    }
}

void file_error( char* filename ) {
    fprintf(stderr, "\n\nUnable to open file %s.\n", filename);
    perror( filename );
    fprintf(stderr, "Leaving RODD.\n");
    fprintf(logfi, "RODD terminated because of error on a file.\n");
    fprintf(logfi, "----------------------------------------------\n");
    fclose(logfi);
    exit(EXIT_FAILURE);
}

void unknown_command( ) {
    fprintf(logfi, "Unknown command found in batch file.\n");
    if (!BATCH) {
	fprintf(stderr, "Unknown command.\n");
    }
    else {
	fprintf(stderr, "\n\nUnknown command found in commmand file.\n");
	fprintf(stderr, "Please fix command file syntax to continue\n.");
	fprintf(stderr, "Terminating RODD.\n");
	fprintf(logfi, "RODD terminated due to unknown command in batch file.\n");
	fprintf(logfi, "----------------------------------------------\n");
	fclose(logfi);
	exit(EXIT_FAILURE);
    }
}

void set_error( ) {
    fprintf(logfi, "Error in set command syntax.\n");
    if (!BATCH) {
	fprintf(stderr, "Unable to set that object.\n");
    }
    else {
	fprintf(stderr, "Error with set command in command file.\n");
	fprintf(stderr, "Please fix command file syntax to continue.\n");
	fprintf(stderr, "Terminating RODD.\n");
	fprintf(logfi, "RODD terminated due to improper set command syntax in batch file.\n");
	fprintf(logfi, "----------------------------------------------\n");
	fclose(logfi);
	exit(EXIT_FAILURE);
    }
}

void addrun_error( char c ) {
    fprintf(logfi, "Error found in addrun syntax.\n");
    switch(c) {
	case ')': fprintf(stderr, "Mismatched parenthesis!\n"); break;
	case '(': fprintf(stderr, "Missing a run number before (!\n"); break;
	case '*': fprintf(stderr, "Improper use of *\n"); break;
	case '-': fprintf(stderr, "Missing a number before dash!\n"); break;
	default: fprintf(stderr, "Illegal character %c.\n", c); 
    }
    if (!BATCH) {
	fprintf(stderr, "Rest of line ignored.\n");
	return;
    }
    else {
	fprintf(stderr, "Error in addrun command in command file\n");
	fprintf(stderr, "Fix command file syntax\n");
	fprintf(stderr, "Terminating RODD.\n");
	fprintf(logfi, "RODD terminated due to incorrect syntax of addrun command.\n");
	fprintf(logfi, "----------------------------------------------\n");
	fclose(logfi);
	exit(EXIT_FAILURE);
    }
}

void size_error( long lower, long upper ) {
    fprintf(logfi, "Start run less than end run, nothing added/\n");
    fprintf(stderr, "The number %l is greater than %l.\n", lower, upper);
    fprintf(stderr, "Check order of numbers when using dash '-' in addrun or delrun\n");
    if (!BATCH) {
	fprintf(stderr, "Command ignored\n");
    }
    else {
	fprintf(stderr, "\n\nSemantic error involving '-' found in command file.\n");
	fprintf(stderr, "Bailing out to allow for a fix.\n");
	fprintf(stderr, "Exiting RODD\n");
	fprintf(logfi, "RODD terminated due to error in batch file.\n");
	fprintf(logfi, "----------------------------------------------\n");
	fclose(logfi);
	exit(EXIT_FAILURE);
    }
}
