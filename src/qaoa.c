/**
 * @author Nicholas Pritchard <21726929@student.uwa.edu.au>
 * @brief The main driving code which actually runs the QAOA
 */
#include "qaoa.h"
#include "ub.h"
#include "uc.h"
#include "state_evolve.h"
#include "reporting.h"
#include "eigen_solve.h"

//TODO Unit test all of the these
/**
 * @brief Checks paramters in the meta-specification for validity
 * @param meta_spec The data-structure containing all relevant fields
 * @warning Number of qubits cannot exceed 31, we are forced to used 32bit integers for intel's FFTW backend
 */
void parameter_checking(qaoa_data_t *meta_spec) {
    //Check machine specification
    if (meta_spec->machine_spec->num_qubits <= 0) {
        fprintf(stderr, "Invalid number of qubits.\n");
        exit(EXIT_FAILURE);
    }
    if (meta_spec->machine_spec->num_qubits > 31) {
        fprintf(stderr, "Too many qubits, numerical stability will fail.\n");
    }
    if (meta_spec->machine_spec->P <= 0) {
        fprintf(stderr, "Invalid amount of decomposition.\n");
        exit(EXIT_FAILURE);
    }

    //Check run specification
    if (meta_spec->run_spec->num_samples <= 0 && meta_spec->run_spec->sampling == true) {
        fprintf(stderr, "Too few samples.\n");
        exit(EXIT_FAILURE);
    }
    if (meta_spec->run_spec->outfile == NULL) {
        fprintf(stderr, "No output location.\n");
        exit(EXIT_FAILURE);
    }

    //Check optimisation specification
    if (meta_spec->opt_spec->max_evals <= 0) {
        fprintf(stderr, "Invalid evaluation count.\n");
        exit(EXIT_FAILURE);
    }

}

/**
 * @brief De-allocates all memory associated with a QAOA simulation
 * @param meta_spec The data-structure containing all relevant fields
 */
void qaoa_teardown(qaoa_data_t *meta_spec){
    mkl_sparse_destroy(meta_spec->ub);
    mkl_free(meta_spec->uc);
    if (!meta_spec->run_spec->restart)
        mkl_free(meta_spec->opt_spec->parameters);
    mkl_free(meta_spec->opt_spec->lower_bounds);
    mkl_free(meta_spec->opt_spec->upper_bounds);
}

//TODO: Include custom optimisation method (not nlopt)
/**
 * @brief Initializes the nlopt optimisation data-structures with provided fields
 * @param meta_spec The data-structure containing all fields
 * @param retain If set, initializer won't allocate or reset parametesr
 * @warning Assumes paramters are of appropriate size and contain relevant values
 */
void optimiser_Initialize(qaoa_data_t *meta_spec, bool retain) {
    int num_params = 2 * meta_spec->machine_spec->P;
    if (meta_spec->run_spec->restricted) {
        num_params++;
    }
    if (!retain) {
        meta_spec->opt_spec->parameters = mkl_calloc((size_t) num_params, sizeof(double), DEF_ALIGNMENT);
        for (int i = 0; i < meta_spec->machine_spec->P; ++i) {
            meta_spec->opt_spec->parameters[i] = (double) PI;
            meta_spec->opt_spec->parameters[i + meta_spec->machine_spec->P] = (double) PI / 2.0;
        }
    } else {
        if (meta_spec->opt_spec->parameters == NULL) {
            fprintf(stderr, "Trying to retain no information\n");
            exit(EXIT_FAILURE);
        }
    }

    meta_spec->opt_spec->lower_bounds = mkl_calloc((size_t) num_params, sizeof(double), DEF_ALIGNMENT);
    meta_spec->opt_spec->upper_bounds = mkl_calloc((size_t) num_params, sizeof(double), DEF_ALIGNMENT);


    for(int i = 0; i < meta_spec->machine_spec->P; ++i){
        meta_spec->opt_spec->upper_bounds[i] = 2 * (double) PI;
        meta_spec->opt_spec->lower_bounds[i] = 0.0;
        meta_spec->opt_spec->upper_bounds[i + meta_spec->machine_spec->P] = (double) PI;
        meta_spec->opt_spec->lower_bounds[i + meta_spec->machine_spec->P] = 0.0;
    }

    meta_spec->opt_spec->optimiser = nlopt_create(meta_spec->opt_spec->nlopt_method, (unsigned int) num_params);

    if (meta_spec->run_spec->restricted) {
        meta_spec->opt_spec->upper_bounds[meta_spec->machine_spec->P] = 2 * (double) PI;
        meta_spec->opt_spec->lower_bounds[meta_spec->machine_spec->P] = 0.0;
        meta_spec->opt_spec->parameters[meta_spec->machine_spec->P] = (double) PI / 2.0;
        nlopt_set_max_objective(meta_spec->opt_spec->optimiser, (nlopt_func) evolve_restricted, (void *) meta_spec);
    } else {
        nlopt_set_max_objective(meta_spec->opt_spec->optimiser, (nlopt_func) evolve, (void *) meta_spec);
    }
    nlopt_set_lower_bounds(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->lower_bounds);
    nlopt_set_upper_bounds(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->upper_bounds);
    nlopt_set_xtol_abs1(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->xtol);
    nlopt_set_ftol_abs(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->ftol);
    nlopt_set_maxeval(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->max_evals);
}

/**
 * @brief The main method which performs a QAOA simulation
 * @param mach_spec Contains the specification of the hypothetial quantum machine
 * @param cost_data Contains information about the cost_function
 * @param opt_spec Contains specification of the classical optimization routine
 * @param run_spec Contains specifcation of the type of simulation to be run
 * @param retain If set, will use parameter values in the opt_spec.
 * @warning If retain set, optimizer will expect values to be pre-initialized
 */
void qaoa(machine_spec_t *mach_spec, cost_data_t *cost_data, optimization_spec_t *opt_spec, run_spec_t *run_spec,
          bool retain) {
    qaoa_data_t meta_spec;
    qaoa_statistics_t statistics;
    MKL_INT ub_nnz;
    statistics.num_evals = 0;
    statistics.best_sample = -INFINITY;
    statistics.best_expectation = -INFINITY;
    meta_spec.qaoa_statistics = &statistics;
    meta_spec.machine_spec = mach_spec;
    meta_spec.run_spec = run_spec;
    meta_spec.opt_spec = opt_spec;
    meta_spec.cost_data = cost_data;

    srand((unsigned) time(0));

    parameter_checking(&meta_spec);

    dsecnd();

    //Initialise UC
    meta_spec.uc = mkl_calloc((size_t) pow(2, meta_spec.machine_spec->num_qubits), sizeof(MKL_Complex16),
                              DEF_ALIGNMENT);
    meta_spec.qaoa_statistics->startTimes[0] = dsecnd();
    meta_spec.qaoa_statistics->startTimes[1] = dsecnd();
    generate_uc(&meta_spec, Cx, mask);
    meta_spec.qaoa_statistics->endTimes[1] = dsecnd();
    if (meta_spec.run_spec->verbose) {
        printf("UC Created\n");
    }
    //Initialise UB
    meta_spec.qaoa_statistics->startTimes[2] = dsecnd();
    ub_nnz = generate_ub(&meta_spec, mask);
    meta_spec.ub_eigenvalue = max_eigen_find(meta_spec.ub);
    //Convert UB to complex values
    convert_ub(&meta_spec, ub_nnz);
    meta_spec.qaoa_statistics->endTimes[2] = dsecnd();
    if (meta_spec.run_spec->verbose) {
        printf("UB Created\n");
    }
    optimiser_Initialize(&meta_spec, retain);
    meta_spec.qaoa_statistics->startTimes[3] = dsecnd();
    meta_spec.qaoa_statistics->term_status = nlopt_optimize(meta_spec.opt_spec->optimiser,
                                                            opt_spec->parameters, &meta_spec.qaoa_statistics->result);
    meta_spec.qaoa_statistics->endTimes[3] = dsecnd();
    meta_spec.qaoa_statistics->endTimes[0] = dsecnd();
    if (meta_spec.run_spec->verbose) {
        printf("Optimisation complete\n");
    }
    //Teardown

    final_report(&meta_spec);
    qaoa_teardown(&meta_spec);
}
