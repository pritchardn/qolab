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

typedef struct {
    double startTimes[4];       // Buffers to hold timing data (total, uc, ub, optimisation)
    double endTimes[4];
    double result;              // The final (previous) function value
    double classical_exp;       // The classical expectation value of random sampling of the cost function
    double random_exp;          // Randomly sampling the entire QAOA domain (may be different to whole state-space)
    double best_result;         // The best cost-value found over the entire optimisation scheme
    int max_value, max_index;   // The maximum value and index in the cost function generated
    int term_status;            // The nlopt termination status
    int num_evals;              // The number of evaluations used by the optimiser
} qaoa_statistics_t;

typedef struct {
    bool timing;        // Do we report timing?
    bool correct;       // Do we report correctness?
    bool report;        // Are we reporting to file?
    bool sampling;      // Are we sampling?
    bool verbose;       // Should we print everything?
    bool restricted;    // Are we running the restricted version of the QAOA (https://arxiv.org/abs/1804.08227)
    int num_samples;    // The number of samples we use
    FILE *outfile;      // The stream we actually write to (can be stdout or a file)
} run_spec_t;

typedef struct {
    int num_qubits;         // The number of qubits in our 'machine'
    int P;                  // The amount of trotterisation
    MKL_INT space_dimension;// The size of the state vector pow(2, qubits)
} machine_spec_t;

typedef struct {
    int nlopt_method;       // The optimisation method used
    int max_evals;          // The maximum number of evaluations permitted
    double xtol;            // The termination criteria for the parameter values
    double ftol;            // The termination criteria for the function output
    double *parameters;     // The parameters to be optimised (2*P length)
    double *lower_bounds;   // The lower bounds for the parameters
    double *upper_bounds;   // The upper bounds for the parameters
    nlopt_opt optimiser;    // The actual optimiser object
    //nlopt_opt local_opt;  //TODO: Support for hybrid multi-optimiser (e.g. MLSL)
} optimisation_spec_t;

/**
 * A meta-structure which contains the information about the entire run
 * This means we can pass a single pointer through many functions but retain access to all elements.
 */
typedef struct{
    machine_spec_t *machine_spec;
    run_spec_t *run_spec;
    cost_data_t *cost_data;
    sparse_matrix_t ub;
    MKL_Complex16 *uc;
    qaoa_statistics_t *qaoa_statistics;
    optimisation_spec_t *opt_spec;
} qaoa_data_t;

#endif