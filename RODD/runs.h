/* runs.h */
#include "ntypes.h"
#include <stdio.h>

/*
 * This is the interface to the list of runs and events 
 * scheduled to be processed when "go" is issued.
 * The set_ routines are called to insert items into the list
 * as they come from the user.
 * Clear_runs and show_listing are also called by the user.
 */

void set_run( int32 );
void set_event( int32, int32 );
void set_multi_events( int32, int32, int32 );
void set_multi_runs( int32, int32 );
void show_listing( FILE* );
void clear_runs( void );
int get_next_node( int32 *, int32 *, int32 * );
void remove_run( int32 );
void multi_remove( int32, int32 );
