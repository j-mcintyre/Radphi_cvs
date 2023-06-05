#ifndef MAP_EDITING_H
#define MAP_EDITING_H

#define MAP_EDITING_STRING_LENGTH 256

#define LOWEST_VALID_RUN 7000
#define HIGHEST_VALID_RUN 9000
#define DETECTOR_NAME_STRING_LENGTH 4
#define MAX_NUMBER_OF_RUNS_BEFORE_BAILOUT 2000

#define UNSIGNED_CHAR_MSB_ON 0x80
#define UNSIGNED_CHAR_ALL_BITS_ON 0xff
#define UNSIGNED_CHAR_N_BITS 8

#define SET_BIT_BIT_ALREADY_SET_RETURN_VALUE 0
#define CLEAR_BIT_NOT_SET_RETURN_VALUE 0

int mark_bit_in_channel_status_map(const char mapfilename[],
				   int detector_number,
				   int channel_number,
				   int first_run,
				   int last_run,
				   unsigned char setbit);
int set_status_bit_in_channel_in_array(unsigned char setbit,
				 int channel_number,
				 unsigned char *array);
int clear_status_bit_in_channel_in_array(unsigned char clearbit,
				 int channel_number,
				 unsigned char *array);
int is_char_bitmask_truly_one_bit(unsigned char bitmask);
int parse_lgdtune_file(FILE *input_file, float *array, int nchannels);

#endif
