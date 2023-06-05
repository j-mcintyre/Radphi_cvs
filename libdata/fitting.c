#include <fitting.h>
#include <string.h>

three_photon_fits_t *create_three_photon_fits_structure(int nfits){
  three_photon_fits_t *output=NULL;
  int i;

  if(nfits<1){
    fprintf(stderr,"create_three_photon_fits_structure called with nfits=%d!!!\n",
	    nfits);
    fprintf(stderr,"Which is less than 1, and is likely an error.\n");
    return(NULL);
  }
  if(nfits>MAX_HYPOTHESIS_FITS){
    fprintf(stderr,"create_three_photon_fits_structure called with nfits=%d!!!\n",
	    nfits);
    fprintf(stderr,"Which is less than the constant MAX_HYPOTHESIS_FITS\n");
    fprintf(stderr,"(defined to be %d in iitypes.h), and is likely an error.\n",
	    MAX_HYPOTHESIS_FITS);
    return(NULL);
  }
  output=(three_photon_fits_t*)malloc(sizeof_three_photon_fits_t(nfits)); 
  if(output==NULL){
    fprintf(stderr,"*create_three_photon_fits_structure:\n");
    fprintf(stderr,"Failed to allocate memory for fits structure, returning NULL\n");
    return(NULL);
  }
  output->nfits=nfits;
  for(i=0;i<nfits;i++){
    output->fit[i].type=FIT_TYPE_UNINITIALIZED;
  }
  return(output);
}

int merge_fit_into_fit_structure(three_photon_fit_t input_fit,
				 three_photon_fits_t *fit_structure,
				 int position,
				 int fit_type){
  if(fit_structure==NULL){
    fprintf(stderr,"merge_fit_into_fit_structure:\n");
    fprintf(stderr,"called with fit_structure set to NULL!\n");
    return(1);
  }
  if(position<0 || position>=fit_structure->nfits){
    fprintf(stderr,"merge_fit_into_fit_structure:\n");
    fprintf(stderr,"called to place fit at position %d\n",position);
    fprintf(stderr,"Which is invalid; given that the number of fits the\n");
    fprintf(stderr,"structure is set up to hold is %d.\n",
	    fit_structure->nfits);
    return(1);
  }
  /* if the input fit is explicitly typed, it's assumed that's the actual
     fit type, and the parameter fit_type is ignored */
  if(fit_type){
    input_fit.type=fit_type;
  }
  else if(input_fit.type==FIT_TYPE_UNINITIALIZED){
    fprintf(stderr,"merge_fit_into_fit_structure:\n");
    fprintf(stderr,"You must either type your fit before merging it,\n");
    fprintf(stderr,"Or pass a type to this function.\n");
    return(1);
  }
  memcpy( (&(fit_structure->fit[position])),(&(input_fit)),
	  sizeof(three_photon_fit_t));
  return(0);
}

void text_dump_fits_structure(three_photon_fits_t *fs){
  int i,j;
  fprintf(stderr,"\ntext_dum_fits_structure:\n");
  fprintf(stderr,"There were %d fits\n",fs->nfits);
  for(i=0;i<fs->nfits;i++){
    fprintf(stderr,"  fit %d:\n",i);
    fprintf(stderr,"    type: %4d\n",fs->fit[i].type);
    fprintf(stderr,"    chi squared: %.3f\n",fs->fit[i].chi_sq);
    fprintf(stderr,"    penalty function: %.3f\n",fs->fit[i].penalty_func);
    fprintf(stderr,"    constraint function: %.3f\n",fs->fit[i].constraint_func);
    fprintf(stderr,"    beam energy: %.2f\n",fs->fit[i].beam_E);
/*
 * Maybe someone should print the appropriate recoil info here [rtj]
 *
    fprintf(stderr,"    recoil type: %d\n",fs->fit[i].recoil.recoil_type);
    fprintf(stderr,"    recoil theta: %.3f\n",fs->fit[i].recoil.theta);
    fprintf(stderr,"    recoil phi:  %.3f\n",fs->fit[i].recoil.phi);
    fprintf(stderr,"    recoil momentum: %.3f\n",fs->fit[i].recoil.momentum_mag);
 */
    for(j=0;j<3;j++){
      fprintf(stderr,"    photon[%d]: E=%5.2f x=%5.2f y=%5.2f z=%5.2f\n",
	      j,
	      fs->fit[i].photon[j].x,
	      fs->fit[i].photon[j].y,fs->fit[i].photon[j].z);
    }
  }
}

int add_selected_fits_to_event(itape_header_t *event,
			       three_photon_fits_t *fs,
			       int N_fits_to_save,
			       int *save_list,
			       int override_full_list_requirement){
  three_photon_fits_t *output_fs=NULL;
  int current_output_position;
  int i,j;

  if(N_fits_to_save<1 || N_fits_to_save>MAX_FITS_TO_SAVE){
    fprintf(stderr,"add_selected_fits_to_event:\n");
    fprintf(stderr,"You have asked to save %d fits, but this\n",
	    N_fits_to_save);
    fprintf(stderr,"is an error!\n");
    return(1);
  }
  if((output_fs=create_three_photon_fits_structure(N_fits_to_save))
     ==NULL){
    fprintf(stderr,"add_selected_fits_to_event:\n");
    fprintf(stderr,"Unable to create local structure for fits!\n");
    return(1);
  }
  current_output_position=0; /* current position of "writing" */
  for(i=0;i<fs->nfits;i++){
    /* i is the fit in fs that we're deciding if we need to write to 
       output_fs */
    for(j=0;j<N_fits_to_save;j++){
      /* j loops through the list of things to save, and decides if 
	 the current fit in the input list (designated by i) matches
	 one of them */
      if(save_list[j]==fs->fit[i].type){
	if(merge_fit_into_fit_structure(fs->fit[i],output_fs,
					current_output_position++,
					fs->fit[i].type)){
	  fprintf(stderr,"add_selected_fits_to_event:\n");
	  fprintf(stderr,"func merge_fit_into_fit_structure\n");
	  fprintf(stderr,"returned an error!\n");
	}
	if(current_output_position>N_fits_to_save){
	  fprintf(stderr,"add_selected_fits_to_event:\n");
	  fprintf(stderr,"DANGER: INTERNAL ERROR!!!\n");
	  fprintf(stderr,"Tried to write %dth fit, even though\n",
		  current_output_position+1);
	  fprintf(stderr,"That %d was the total number of output fits!\n",
		  current_output_position);
	  return(1);
	}
      }
    }
  }
  if(!override_full_list_requirement)
    if(current_output_position!=N_fits_to_save){
      fprintf(stderr,"WARNING!!!\n");
      fprintf(stderr,"add_selected_fits_to_event:\n");
      fprintf(stderr,"Went all the way through the input list of fits,\n");
      fprintf(stderr,"but did not find the specified number to save.\n");
    }
  if(add_fits_group_to_event(event,output_fs)){
    fprintf(stderr,"add_selected_fits_to_event:\n");
    fprintf(stderr,"output list of fits sorted without problems,\n");
    fprintf(stderr,"but add_fits_group_to_event returned an error.\n");
    return(1);
  }
  free(output_fs);
  return(0);
}
			       
int add_fits_group_to_event(itape_header_t *event,
			    three_photon_fits_t *fs){
  void *return_pointer;
  size_t group_size;
  group_size=sizeof_three_photon_fits_t(fs->nfits);
  return_pointer=(data_addGroup(event,BUFSIZE,GROUP_3_PHOTON_FITS,0,
				group_size));
  if(return_pointer==NULL){
    fprintf(stderr,"add_fits_group_to_event:\n");
    fprintf(stderr,"In attempting to add the fits group to the event,\n");
    fprintf(stderr,"data_addGroup has returned NULL.  I expect that's probably\n");
    fprintf(stderr,"an error condition.\n");
    return(1);
  }
  memcpy(return_pointer,fs,group_size);
  return(0);
}

int write_3photon_dst(itape_header_t *event,
			  FILE *output_file){
  /* this function "saves" groups (decides which ones need to be kept)
     does a "clean" (removes the groups not to be kept), and then
     writes the data to a file */
#define N_GROUPS_TO_KEEP 4
  int return_value;
  int n_groups_to_keep=N_GROUPS_TO_KEEP;
  int keeper_array_groups[N_GROUPS_TO_KEEP]={GROUP_3_PHOTON_FITS,
					     GROUP_LGD_CLUSTERS,
					     GROUP_PHOTONS,
                                             GROUP_RAW};
  int keeper_array_subgroups[N_GROUPS_TO_KEEP]={0,0,0};
  if(return_value=data_saveGroups(event,n_groups_to_keep,
			       keeper_array_groups,
				  keeper_array_subgroups)){
    fprintf(stderr,"write_out_3photon_dst:\n");
    fprintf(stderr,"return value from data_saveGroups was %d\n",
	    return_value);
  }
  if(data_clean(event)){
    fprintf(stderr,"write_out_3photon_dst:\n");
    fprintf(stderr,"data_clean returned non-zero.\n");
    return(1);
  }
  if((return_value=data_write(fileno(output_file),event))!=DATAIO_OK){
    fprintf(stderr,"write_out_3photon_dst:\n");
    fprintf(stderr,"data_write returned %d.\n",return_value);
    fprintf(stderr,"which is not equal to DATAIO_OK.\n");
    return(1);
  }
  return(0);
}
