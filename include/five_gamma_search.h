#ifndef FIVE_GAMMA_SEARCH_H
#define FIVE_GAMMA_SEARCH_H

/* this file contains utility structures that may be useful in doing
   searches with combinations of 5 photons.  */

#define FIVE_GAMMA_N_GROUPINGS 15
#define FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS 30
#define FIVE_GAMMA_PAIRCENTERED_PAIR_INTERVAL 3
#define FIVE_GAMMA_PAIRCENTERED_N_PAIRS 10

/* external versions */
extern int g5_grp_bat[FIVE_GAMMA_N_GROUPINGS];
extern int g5_grp_pr0_g0[FIVE_GAMMA_N_GROUPINGS];
extern int g5_grp_pr0_g1[FIVE_GAMMA_N_GROUPINGS];
extern int g5_grp_pr1_g0[FIVE_GAMMA_N_GROUPINGS];
extern int g5_grp_pr1_g1[FIVE_GAMMA_N_GROUPINGS];
int g5_pr_pr0_g0[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS];
int g5_pr_pr0_g1[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS];
int g5_pr_pr1_g0[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS];
int g5_pr_pr1_g1[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS];
int g5_pr_bat[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS];


/* real definitions */
#ifdef FILE_CONTAINING_MAIN

int g5_grp_bat[FIVE_GAMMA_N_GROUPINGS]=
{0,0,0,1,1,1,2,2,2,3,3,3,4,4,4};

int g5_grp_pr0_g0[FIVE_GAMMA_N_GROUPINGS]=
{1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};

int g5_grp_pr0_g1[FIVE_GAMMA_N_GROUPINGS]=
{2,3,4,2,3,4,1,3,4,1,2,4,1,2,3};

int g5_grp_pr1_g0[FIVE_GAMMA_N_GROUPINGS]=
{3,2,2,3,2,2,3,1,1,2,1,1,2,1,1};

int g5_grp_pr1_g1[FIVE_GAMMA_N_GROUPINGS]=
{4,4,3,4,4,3,4,4,3,4,4,2,3,3,2};

int g5_pr_pr0_g0[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS]=
{0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 1,1,1, 1,1,1, 2,2,2, 2,2,2, 3,3,3};

int g5_pr_pr0_g1[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS]=
{1,1,1, 2,2,2, 3,3,3, 4,4,4, 2,2,2, 3,3,3, 4,4,4, 3,3,3, 4,4,4, 4,4,4};

int g5_pr_pr1_g0[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS]=
{2,2,3, 1,1,3, 1,1,2, 1,1,2, 0,0,4, 0,0,4, 0,0,3, 0,0,4, 0,0,3, 0,0,2};

int g5_pr_pr1_g1[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS]=
{3,4,4, 3,4,4, 2,4,4, 2,3,3, 4,3,3, 4,2,2, 3,2,2, 4,1,1, 3,1,1, 2,1,1};

int g5_pr_bat[FIVE_GAMMA_N_PAIRCENTERED_GROUPINGS]=
{4,3,2, 4,3,1, 4,2,1, 3,2,1, 3,4,0, 2,4,0, 2,3,0, 1,4,0, 1,3,0, 1,2,0};


#endif

#endif
