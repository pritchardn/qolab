#ifndef GRAPHSIMILARITY_GLOBALS_H
#define GRAPHSIMILARITY_GLOBALS_H


#include <stdio.h>
#include <mkl.h>
#include <stdbool.h>
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
} run_spec_t;

typedef struct {
    int num_qubits;
    int P;
    MKL_INT space_dimension;
} machine_spec_t;

typedef struct {
    int nlopt_method;
    double xtol;
    double ftol;
} optimisation_spec_t;


typedef struct{
    machine_spec_t *machine_spec;
    run_spec_t *run_spec;
    sparse_matrix_t ub;
    MKL_Complex16 *uc;
    cost_data_t *cost_data;
    qaoa_statistics_t *qaoa_statistics;
}qaoa_data_t;

#endif