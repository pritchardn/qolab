/**
 * @author Nicholas Pritchard <21726929@student.uwa.edu.au>
 * The main driving code which actually runs the QAOA
 */
#include "qaoa.h"
#include "ub.h"
#include "uc.h"
#include "state_evolve.h"
#include "reporting.h"

//TODO Unit test all of the these
void parameter_checking(qaoa_data_t *meta_spec) {
    //Check machine specification
    if (meta_spec->machine_spec->num_qubits <= 0) {
        fprintf(stderr, "Invalid number of qubits.\n");
        exit(EXIT_FAILURE);
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

void qaoa_teardown(qaoa_data_t *meta_spec){
    mkl_sparse_destroy(meta_spec->ub);
    mkl_free(meta_spec->uc);
    mkl_free(meta_spec->opt_spec->parameters);
    mkl_free(meta_spec->opt_spec->lower_bounds);
    mkl_free(meta_spec->opt_spec->upper_bounds);
}

//TODO: Include custom optimisation method (not nlopt)
void optimiser_initialise(qaoa_data_t *meta_spec){
    meta_spec->opt_spec->parameters = mkl_calloc((size_t)2 * meta_spec->machine_spec->P, sizeof(double), DEF_ALIGNMENT);
    meta_spec->opt_spec->lower_bounds = mkl_calloc((size_t)2 * meta_spec->machine_spec->P, sizeof(double), DEF_ALIGNMENT);
    meta_spec->opt_spec->upper_bounds = mkl_calloc((size_t)2 * meta_spec->machine_spec->P, sizeof(double), DEF_ALIGNMENT);

    for(int i = 0; i < meta_spec->machine_spec->P; ++i){
        meta_spec->opt_spec->upper_bounds[i] = 2 * (double)M_PI;
        meta_spec->opt_spec->lower_bounds[i] = 0.0;
        meta_spec->opt_spec->upper_bounds[i+meta_spec->machine_spec->P] = (double)M_PI;
        meta_spec->opt_spec->lower_bounds[i+meta_spec->machine_spec->P] = 0.0;
        meta_spec->opt_spec->parameters[i] = (double)M_PI;
        meta_spec->opt_spec->parameters[i + meta_spec->machine_spec->P] = (double)M_PI/2.0;
    }

    meta_spec->opt_spec->optimiser = nlopt_create(meta_spec->opt_spec->nlopt_method, (unsigned int)2*meta_spec->machine_spec->P);
    nlopt_set_max_objective(meta_spec->opt_spec->optimiser, (nlopt_func)evolve, (void*)meta_spec);
    nlopt_set_lower_bounds(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->lower_bounds);
    nlopt_set_upper_bounds(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->upper_bounds);
    nlopt_set_xtol_abs1(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->xtol);
    nlopt_set_ftol_abs(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->ftol);
    nlopt_set_maxeval(meta_spec->opt_spec->optimiser, meta_spec->opt_spec->max_evals);
}

void qaoa(machine_spec_t *mach_spec, cost_data_t *cost_data, optimisation_spec_t *opt_spec, run_spec_t *run_spec){
    qaoa_data_t meta_spec;
    qaoa_statistics_t statistics;
    statistics.num_evals = 0;
    statistics.best_result = -INFINITY;
    meta_spec.qaoa_statistics = &statistics;
    meta_spec.machine_spec = mach_spec;
    meta_spec.run_spec = run_spec;
    meta_spec.opt_spec = opt_spec;
    meta_spec.cost_data = cost_data;

    srand((unsigned) time(0));

    parameter_checking(&meta_spec);

    dsecnd();

    //Initialise UC
    meta_spec.uc = mkl_malloc((size_t)pow(2, meta_spec.machine_spec->num_qubits) * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    meta_spec.qaoa_statistics->startTimes[0] = dsecnd();
    meta_spec.qaoa_statistics->startTimes[1] = dsecnd();
    generate_uc(&meta_spec, Cx);
    meta_spec.qaoa_statistics->endTimes[1] = dsecnd();
    if (meta_spec.run_spec->verbose) {
        printf("UC Created\n");
    }
    //Initialise UB
    meta_spec.qaoa_statistics->startTimes[2] = dsecnd();
    generate_ub(&meta_spec);
    meta_spec.qaoa_statistics->endTimes[2] = dsecnd();
    if (meta_spec.run_spec->verbose) {
        printf("UB Created\n");
    }
    optimiser_initialise(&meta_spec);
    meta_spec.qaoa_statistics->startTimes[3] = dsecnd();
    meta_spec.qaoa_statistics->term_status = nlopt_optimize(meta_spec.opt_spec->optimiser, opt_spec->parameters, &meta_spec.qaoa_statistics->result);
    meta_spec.qaoa_statistics->endTimes[3] = dsecnd();
    meta_spec.qaoa_statistics->endTimes[0] = dsecnd();
    if (meta_spec.run_spec->verbose) {
        printf("Optimisation complete\n");
    }
    //Teardown

    final_report(&meta_spec);

    qaoa_teardown(&meta_spec);
}
