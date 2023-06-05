/* runs.c */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "runs.h"
#include "errors.h"
#include "ntypes.h"
#include "global.h"


typedef struct to_process {
  int32 the_run;             /* the run number to process */
  int32 first_event;        /* the first event in this set to process inclusive */
  int32 last_event;         /* the last event to process for this run inclusive */
  struct to_process *next; /* pointer to next node */
} node;

static node *head = NULL;	    /* head of the list of run nodes */

/* build_node does memory allocation and returns a pointer to the
 * new node containing the data that gets passed in
 */

static
node* build_node( int32 run, int32 low, int32 up ) {
    node* temp;
    if ((temp = (node *)malloc(sizeof(node))) == NULL)
      malloc_error(__FILE__, __LINE__);
    temp->the_run = run;
    temp->first_event = low;
    temp->last_event = up;
    temp->next = NULL;
    return(temp);
}

/* report informs the user what got inserted into the list */

static
void report( node* inserted ) {
    if ( inserted->first_event == inserted->last_event )
	fprintf(logfi, "Event number %d added for run %i.\n", 
	    inserted->first_event,  inserted->the_run);
    else if ( (inserted->first_event == 1) && (inserted->last_event == INT_MAX))
	fprintf(logfi, "All events of run %i added.\n", inserted->the_run);
    else
	fprintf(logfi, "Events %d through %d added for run %i.\n", 
	    inserted->first_event, inserted->last_event, inserted->the_run);
    if(!BATCH) {
	if ( inserted->first_event == inserted->last_event )
	    fprintf(stdout, "Event number %d added for run %i.\n", 
		inserted->first_event,  inserted->the_run);
	else if ( (inserted->first_event == 1) && (inserted->last_event == INT_MAX))
	    fprintf(stdout, "All events of run %i added.\n", inserted->the_run);
	else
	    fprintf(stdout, "Events %d through %d added for run %i.\n", 
		inserted->first_event, inserted->last_event, inserted->the_run);
    }
}

/* insert does a sorted insert into the list of runs.
 * the new data node is sorted first by run number.
 * if the new data node's run number is already on this list, we
 * then sort according to the starting event number
 */

static void 
insert( node* data ) {
  node* prev = NULL;
  node* curr = head;
  
  if (head == NULL) head = data;		    /*empty list*/
  else {
    while ((curr != NULL) && (curr->the_run < data->the_run)) {
      prev = curr;
      curr = curr->next;
    } /* while */

    if (curr == NULL) prev->next = data;	    /*end of list*/
    else if (curr->the_run != data->the_run) {      /*new run #*/
      if (curr == head) {			    /*at head*/
	data->next = head;
	head = data;
      } /*if curr == head */
      else {				            /*elsewhere*/
	prev->next = data;
	data->next = curr;
      }
    } /* else if */
    else {				/*same run #, sort by event*/
      while((curr != NULL) && (curr->the_run == data->the_run)) {
	if ((data->first_event < curr->first_event) && (data->last_event < curr->first_event)) {
	  /* insert here! */
	  if (curr == head) { /* at head */
	    data->next = head;
	    head = data;
	    data = NULL;
	  }
	  else {
	    prev->next = data;
	    data->next = curr;
	    data = NULL;
	  }
	  break;
	}
	else if (data->first_event > curr->last_event) {
	  prev = curr; curr = curr->next;
	}
	else if ((data->first_event < curr->first_event) && (data->last_event <= curr->last_event)) {
	  curr->first_event = data->first_event;
	  report(data); free(data); return;
	  break;
	}
	else if ((data->first_event >= curr->first_event) && (data->last_event <= curr->last_event)) {
	  report(data); free(data); return;
	  break;
	}
	else if ((data->first_event >= curr->first_event) && (data->last_event > curr->last_event)) {
	  curr->last_event = data->last_event;
	  report(data); free(data); return;
	  break;
	}
	else if ((data->first_event <= curr->first_event) && (data->last_event >= curr->last_event)) {
	  curr->first_event = data->first_event;
	  curr->last_event = data->last_event;
	  report(data); free(data); return;
	  break;
	}
	else fprintf(stdout, "A GOOF!\n");
      } /* end while run#s equal */
      if (data != NULL) {  /* still need to insert! */
	if (curr == head) {
	  data->next = head;
	  head = data;
	}
	else {
	  prev->next = data;
	  data->next = curr;
	}
      }
    } /* else same run #, event sort */
  } /* else, ie not at head of list */
  report(data);
} /* insert */

static void 
delete_node( int32 arun ) {
  node* curr = head;
  node* prev = NULL;
  int count = 0;
  
  while(curr != NULL) {
    while( (curr != NULL) && (curr->the_run != arun ) ) {
      prev = curr;
      curr = curr->next;
    }
    if (curr == NULL) break;	    /* done with list */
    if (curr == head) {		    /* remove the head */
      head = head->next;
      free(curr);
      curr = head;
    }
    else {			    /* middle of list */
      prev->next = curr->next;
      curr->next = NULL;
      free(curr);
      curr = prev->next;
    }
    count++;
  }
  if (count > 0) 
    fprintf(logfi, "\n%i sets of run number %i data removed", count, arun);
    if(!BATCH) fprintf(stdout, "\n%i sets of run number %i data removed", count, arun);
  else
    if(!BATCH) fprintf(stdout, "\nRun number %i not found in list, no data changed.", arun);
    fprintf(logfi, "\nRun number %i not found in list, no data changed.", arun);
}

/* PUBLIC FUNCTIONS */
    
void set_run( int32 run ) {
    insert(build_node(run, 1, INT_MAX));
}

void set_event( int32 run, int32 event ) {
    insert(build_node(run, event, event));
}

void set_multi_events( int32 run, int32 lowevent, int32 topevent ) {
    insert(build_node(run, lowevent, topevent));
}

void set_multi_runs( int32 lowrun, int32 toprun ) {
    int i;
    for( i = lowrun; i <= toprun; i++) 
	insert(build_node(i, 1, INT_MAX));  
}

void show_listing(FILE* fp) {
  node* curr = head;
  if (curr == NULL) fprintf(fp, "\nNo runs set to be processed");
  while ( curr != NULL ) {
    fprintf(fp, "Run %i. ", curr->the_run);
    if (curr->first_event == curr->last_event)
	fprintf(fp, "Event %i.\n", curr->first_event);
    else if (curr->last_event == INT_MAX)
	fprintf(fp, "Events %i through end of run.\n",  curr->first_event);
    else fprintf(fp, "Events %i through %i.\n",  curr->first_event, curr->last_event);
    curr = curr->next;
  }
}

void clear_runs() {
  node* curr = head;
  node* to_free = NULL;
  while (curr != NULL) {
    to_free = curr;
    curr = curr->next;
    free(to_free);
  }
  head = NULL;
  fprintf(logfi, "List cleared\n");
  if(!BATCH) fprintf(stdout, "List cleared\n");
}

/*
 * remove_run removes all instances of the run from the list
 */
 
void remove_run( int32 run ) {
    delete_node(run);
}

void multi_remove( int32 lower, int32 upper ) {
    int i;
    for (i = lower; i <= upper; i++)
	delete_node(i);
}

int get_next_node( int32 *run, int32 *lower, int32 *upper) {
    if (head != NULL) {
	node* temp = head->next;
	*run = head->the_run;
	*lower = head->first_event;
	*upper = head->last_event;
	free(head);
	head = temp;
	return 1;
    }
    else return 0;
}

