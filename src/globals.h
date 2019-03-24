/**
 * @author Nicholas Pritchard
 * @date 28/06/2018
 */
#ifndef GRAPHSIMILARITY_GLOBALS_H
#define GRAPHSIMILARITY_GLOBALS_H

#include <stdio.h>
#include <mkl.h>
#include <stdbool.h>
#include <nlopt.h>
#include "problem_code.h"

#define DEF_ALIGNMENT 64
#define PI 3.1415926535


//Mathematical
MKL_INT factorial(int n);

//Administrative
void mkl_error_parse(int error, FILE *stream);
void check_alloc(void *pointer);

void move_params(int P, double *parameters);

void move_params_restricted(int P, double *parameters);

/*! Contains run-time statistics */
typedef struct {
    double startTimes[4];       /**< Buffers to hold timing data (total, uc, ub, optimisation) */
    double endTimes[4];
    double result;              /**< The final (previous) function value */
    double classical_exp;       /**< The classical expectation value of random sampling of the cost function */
    double random_exp;          /**< Randomly sampling the entire QAOA domain (may be different to whole state-space) */
    double best_sample;    /**< The best expectation value found */
    double best_expectation;
    int max_value, max_index;   /**< The maximum value and index in the cost function generated */
    nlopt_result term_status;   /**< The nlopt termination status */
    int num_evals;              /**< The number of evaluations used by the optimiser */
} qaoa_statistics_t;

/*! Defines run-time parameters on what to report and the type of algorithm simulated */
typedef struct {
    bool timing;        /**< Do we report timing? */
    bool correct;       /**< Do we report correctness? */
    bool report;        /**< Are we reporting to file? */
    bool sampling;      /**< Are we sampling? */
    bool verbose;       /**< Should we print everything? */
    bool restricted;    /**< Are we running the restricted version of the QAOA? (https://arxiv.org/abs/1804.08227) */
    bool restart;       /**< If set, the simulation will retain parameter information between calls to the simulation */
    int num_samples;    /**< The number of samples we use */
    FILE *outfile;      /**< The stream we actually write to (can be stdout or a file) */
} run_spec_t;

/*! Defines the size of the quantum machine to be simulated */
typedef struct {
    int num_qubits;          /**< The number of qubits in our 'machine' */
    int P;                   /**< The amount of trotterisation */
    MKL_INT space_dimension; /**< The size of the state vector pow(2, qubits) */
} machine_spec_t;

/*! Specifies the classical optimisation scheme */
typedef struct {
    int nlopt_method;       /**< The optimisation method used */
    int max_evals;          /**< The maximum number of evaluations permitted */
    double xtol;            /**< The termination criteria for the parameter values */
    double ftol;            /**< The termination criteria for the function output */
    double *parameters;     /**< The parameters to be optimised (2*P length) */
    double *lower_bounds;   /**< The lower bounds for the parameters */
    double *upper_bounds;   /**< The upper bounds for the parameters */
    nlopt_opt optimiser;    /**< The actual optimiser object */
    //nlopt_opt local_opt;  //TODO: Support for hybrid multi-optimiser (e.g. MLSL)
} optimization_spec_t;

/*! A meta-structure which contains the information about the entire run
 *
 * This means we can pass a single pointer through many functions but retain access to all elements. */
typedef struct{
    machine_spec_t *machine_spec;        /**< Defines the quantum machine */
    run_spec_t *run_spec;               /**< Defines run-time parameters */
    cost_data_t *cost_data;             /**< Contains the problem-dependent information */
    sparse_matrix_t ub;                 /**< The data structure containing the driver hamiltonian */
    double ub_eigenvalue;               /**< The leading eigenvalue of the UB matrix */
    MKL_Complex16 *uc;                  /**< The data structure containing the cost function */
    qaoa_statistics_t *qaoa_statistics; /**< Contains run-time statistics */
    optimization_spec_t *opt_spec;      /**< Specifies the classical optimisation scheme */
} qaoa_data_t;

#endif