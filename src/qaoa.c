#include "qaoa.h"
#include "ub.h"
#include "uc.h"
#include "problem_code.h"
#include "state_evolve.h"
#include <stdio.h>
#include <nlopt.h>


void qaoa_teardown(qaoa_data_t *meta_spec){
    mkl_sparse_destroy(meta_spec->ub);
    mkl_free(meta_spec->uc);
}

void optimiser_initialise(qaoa_data_t *meta_spec){
    meta_spec->opt_spec->parameters = mkl_calloc((size_t)2*meta_spec->machine_spec->P, sizeof(double), DEF_ALIGNMENT);
    meta_spec->opt_spec->lower_bounds = mkl_calloc((size_t)2*meta_spec->machine_spec->P, sizeof(double), DEF_ALIGNMENT);
    meta_spec->opt_spec->upper_bounds = mkl_calloc((size_t)2*meta_spec->machine_spec->P, sizeof(double), DEF_ALIGNMENT);

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
    meta_spec.qaoa_statistics = &statistics;
    meta_spec.machine_spec = mach_spec;
    meta_spec.run_spec = run_spec;
    meta_spec.opt_spec = opt_spec;


    //Initialise UC
    meta_spec.uc = mkl_malloc((size_t)pow(2, meta_spec.machine_spec->num_qubits) * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    generate_uc(&meta_spec, cost_data, Cx);
    //Initialise UB
    generate_ub(&meta_spec);
    //mkl_sparse_print(&meta_spec.ub, stdout);
    //Trial feval
    optimiser_initialise(&meta_spec);
    double result=0.0;
    nlopt_optimize(meta_spec.opt_spec->optimiser, opt_spec->parameters, &result);
    printf("%f\n", result);

    //Teardown
    mkl_free(opt_spec->parameters);
    mkl_free(opt_spec->lower_bounds);
    mkl_free(opt_spec->upper_bounds);
    qaoa_teardown(&meta_spec);
}
