#ifndef TRIGGERCONTROL_H
#define TRIGGERCONTROL_H

#define TRIGGER_INFO_TYPE_NONE 0
#define TRIGGER_INFO_TYPE_TF   1
#define TRIGGER_INFO_TYPE_INT  2
#define MAX_LINE_LENGTH 256


typedef struct{
  int type;
  int in_latch;
  unsigned long value;
  int bypassed;
}trig_info_t;

#define N_TRIGGER_SETTINGS 7
typedef struct{
  int set_up;
  int latch_slot;
  unsigned long latch_mask;
  int mlu_slot;
  trig_info_t BSD_R;
  trig_info_t BSD_L;
  trig_info_t BSD_S;
  trig_info_t TAG;
  trig_info_t UPV;
  trig_info_t L2;
  trig_info_t L3;
}master_trigger_t;

void trigger_utilities_global_init(master_trigger_t **input_trigger);
int get_crate_slot_of_module(char *infile_name,
			     char *mod_name,
			     int *crate,
			     int *slot);
int find_module_colon(char *line, int *position);
int find_word_colon(char *line, char *word, int length, int *position);
int parse_out_crate_slot(char *line, int *crate, int *slot);
int parse_dotted_int(char *line, int *num);
int clip_end_by_one(char *string);
int is_a_bit(unsigned long val);
int strip_out_value_from_line(char *line,
			      char *indicator, unsigned long *data);
int strip_out_logical_from_line(char *line, char *indicator, 
				int *data);
int parse_xrdb_for_info(char *xrdb_file_name,
			master_trigger_t *trigger_param);
int is_line_commented(char *line);
int print_trigger_factoids(FILE *outfile, master_trigger_t *trigger);
int trig_print_decimal_line(FILE *outfile, char *line, int position, 
			    unsigned long value);	 
int trig_print_hex_line(FILE *outfile, char *line, int position, 
			    unsigned long value);	 
int make_new_camac_file(char *old_file_name,
			char *new_file_name,
			master_trigger_t *spec);
void re_synch_trigger(master_trigger_t *trig);
int get_logical_from_console(const char *text, int *logical);
int update_trigger(master_trigger_t *trigger, int *trigger_modified);
int modify_trigger(master_trigger_t *trigger);
#endif
