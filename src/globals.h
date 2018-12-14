#ifndef GRAPHSIMILARITY_GLOBALS_H
#define GRAPHSIMILARITY_GLOBALS_H


#include <stdio.h>
#include <mkl.h>
#include <stdbool.h>
#include <nlopt.h>
#include "problem_code.h"

#define DEF_ALIGNMENT 64


//Mathematical
MKL_INT factorial(int n);

//Administrative
void mkl_error_parse(int error, FILE *stream);
void check_alloc(void *pointer);

//Debugging
void mkl_sparse_print(sparse_matrix_t *matrix, FILE *stream);

typedef struct {
    double startTimes[5];
    double endTimes[5];
    double result, classical_exp, random_exp;
    int max_value, max_index;
    int term_status;
} qaoa_statistics_t;

typedef struct {
    bool timing;
    bool correct;
    bool report;
    bool sampling;
    int num_samples;
} run_spec_t;

typedef struct {
    int num_qubits;
    int P;
    MKL_INT space_dimension;
} machine_spec_t;

typedef struct {
    int nlopt_method;
    int max_evals;
    double xtol;
    double ftol;
    double *parameters;
    double *lower_bounds;
    double *upper_bounds;
    nlopt_opt optimiser;
    //nlopt_opt local_opt;
} optimisation_spec_t;


typedef struct{
    machine_spec_t *machine_spec;
    run_spec_t *run_spec;
    sparse_matrix_t ub;
    MKL_Complex16 *uc;
    qaoa_statistics_t *qaoa_statistics;
    optimisation_spec_t *opt_spec;
} qaoa_data_t;

#endif