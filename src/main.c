/**
 * @author Nicholas Pritchard <21726929@student.uwa.edu.au>
 * @version 1.0
 */

#include "qaoa.h"
#include "graph_utils.h"
#include <mathimf.h>

void init_graph(int num_qubits, cost_data_t *cost_data) {
    int i, j;
    cost_data->graph = mkl_calloc((size_t) num_qubits * num_qubits, sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(cost_data->graph);
    for (i = 0; i < num_qubits - 1; ++i) {
        j = i + 1;
        cost_data->graph[i * num_qubits + j] = 1;
        cost_data->graph[j * num_qubits + i] = 1;
    }
    cost_data->graph[num_qubits - 1] = 1;
    cost_data->graph[(num_qubits - 1) * num_qubits] = 1;
    print_graph(cost_data, stdout);
}

/*
 * An example main file which runs our example solution
 * All parameters which must be specified are in this example.
 */
int main(int argc, char *argv[]){
    run_spec_t run_spec;
    run_spec.correct = true;
    run_spec.report = false;
    run_spec.timing = true;
    run_spec.sampling = false;
    run_spec.verbose = false;
    run_spec.restricted = true;
    run_spec.num_samples = 100;
    run_spec.outfile = stdout;

    machine_spec_t mach_spec;
    mach_spec.num_qubits = 4;
    mach_spec.P = 1;
    mach_spec.space_dimension = (MKL_INT)pow(2, mach_spec.num_qubits);

    optimisation_spec_t opt_spec;
    opt_spec.ftol = 1e-16;
    opt_spec.xtol = 1e-16;
    opt_spec.nlopt_method = NLOPT_LN_NELDERMEAD;
    opt_spec.max_evals = 200 * mach_spec.P;

    cost_data_t cost_data;
    cost_data.cx_range = mach_spec.space_dimension;
    cost_data.x_range = mach_spec.space_dimension;
    cost_data.num_vertices = mach_spec.num_qubits;
    cost_data.graph = mkl_malloc(sizeof(MKL_INT) * mach_spec.num_qubits * mach_spec.num_qubits, DEF_ALIGNMENT);
    init_graph(mach_spec.num_qubits, &cost_data);
    qaoa(&mach_spec, &cost_data, &opt_spec, &run_spec);
    mkl_free(cost_data.graph);
    return 0;
}
