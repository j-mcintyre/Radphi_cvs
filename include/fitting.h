#include<stdio.h>
#include<stdarg.h>
#include<disData.h>
#include<itypes.h>
#include<iitypes.h>
#include<makeHits.h>
#include<dataIO.h>

#ifndef FITTING_H
#define FITTING_H

#define MAX_FITS_TO_SAVE 100

three_photon_fits_t *create_three_photon_fits_structure(int nfits);
int merge_fit_into_fit_structure(three_photon_fit_t input_fit,
				 three_photon_fits_t *fit_structure,
				 int position,
				 int fit_type);
int add_fits_group_to_event(itape_header_t *event,
			    three_photon_fits_t *fs);
int writeadd_selected_fits_to_event_3photon_dst(itape_header_t *event,
		  FILE *output_file);
int do_fit_stuff(void);
void text_dump_fits_structure(three_photon_fits_t *fs);
int add_selected_fits_to_event(itape_header_t *event,
			       three_photon_fits_t *fs,
			       int N_fits_to_save,
			       int *save_list,
			       int override_full_list_requirement);


#endif
