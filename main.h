#ifndef _EDF_MAIN_H
#define _EDF_MAIN_H

#include "edf_hull.h"
#include "ts_lib.h"

#include <argp.h>
#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LEN(v) (sizeof(v) / sizeof(v[0]))

#define MAX_SIZE_STR_FILENAME 500
#define LENGTH_COMMAND 500
#define LEN_TASK_GENERATING_DL 100
#define NO_SEED -1

#define EDFH_SEED 140
#define EDFH_NUM_TASKS 141
#define EDFH_PERIOD_MIN 142
#define EDFH_PERIOD_MAX 143
#define EDFH_PHASING 144
#define EDFH_RELATIVE_DL_AVG 145
#define EDFH_RELATIVE_DL_VAR 146
#define EDFH_EPS 147
#define EDFH_N_REPEAT 148
#define EDFH_DL_MIN 149
#define EDFH_DL_MAX 150
#define EDFH_SEED_POINT 151

/*
 * [FILE MODE] Option 'i': apply the convex hull point reduction to the task set
 * read from stdin (from text file via stdin redirection). If the
 * option is "iv", then verbose output
 */
void main_file_mode();

/*
 * [RANDOM MODE] Option 's': apply the convex hull point reduction to a single
 * task set randomly generated by the random seed and the other settings passed
 * via stdin . If the option is "cv", then verbose output.
 * If combined with option "num-repeat", then the convex hull point reduction
 * will be applied to num_repeat task sets randomly generated by the random seed
 * and the other.
 */
void main_rand_mode();

/* [RANDOM MODE] Option 'e': apply the convex hull point reduction to one task
 * set randomly generated by the random seed and the other settings passed via
 * stdin. Plus*/
void main_constraints_info();

/**
 * [RANDOM MODE] Option 'num-repeat' is included: apply the convex hull point
 * reduction to multiple task sets and print results on csv file in directory
 * '../datasets'
 */
void iterate_random_mode();

/**
 * Create a descriptive filename based on the random settings
 */
void create_descriptive_filename(char *filename);

void print_rand_setup(ts_rand_t rand_setup);

int verify_arguments_random_mode();

/**
 * Print the results of the procedure on the terminal
 */
void edf_print_stats(edf_points_t *my_points, const ts_t *my_task_set,
                     double time_points, double time_qhull);

/**
 * Print the results of the procedure on csv file positioned in the directory
 * ../datasets/new/
 */
void edf_print_stats_on_csv(const ts_rand_t *settings, edf_points_t *my_points,
                            const ts_t *my_task_set);

/**
 * Print the additional information of the
 */

void edf_print_additional_info_on_csv(edf_points_t *my_points,
                                      const ts_t *my_task_set, int seed);

/**
 * Verify if the argument is a valid integer
 */
int verify_arg_int(char *arg);

/**
 * This function is used to get the command line arguments in order to print it
 * on the file CSV containing the results of the execution of the program.
 * The purpose of this is to be able to reproduce the results of the execution
 * of the program on a randomly generated task set.
 */
void get_command(int argc, char **argv);

#endif