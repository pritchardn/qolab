/**
 * @author Nicholas Pritchard <21726929@student.uwa.edu.au>
 * @version 1.0
 */

#include "qaoa.h"
#include "graph_utils.h"
#include "eigen_solve.h"
#include <mathimf.h>

/*
 * An example main file which runs our example solution
 * All parameters which must be specified are in this example.
 */
int main(int argc, char *argv[]){
    MKL_INT ia[5] = {1, 3, 5, 7, 9};                         /* ia array from CSR format */
    MKL_INT ja[8] = {1, 2, 1, 2, 3, 4, 3, 4};                   /* ja array from CSR format */
    double a[8] = {6.0, 2.0, 2.0, 3.0, 2.0, -1.0, -1.0, 2.0}; /* val array from CSR format */
    sparse_matrix_t A;
    sparse_status_t status;
    status = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ONE, 4, 4, ia, ia + 1, ja, a);
    if (status != 0) {
        printf("error building matrix\n");
        exit(EXIT_FAILURE);
    }
    max_eigen_find(A);
    exit(EXIT_SUCCESS);
    run_spec_t run_spec;
    run_spec.correct = true;
    run_spec.report = false;
    run_spec.timing = true;
    run_spec.sampling = true;
    run_spec.verbose = false;
    run_spec.restricted = true;
    run_spec.num_samples = 100;
    run_spec.outfile = stdout;

    machine_spec_t mach_spec;
    mach_spec.num_qubits = 4;
    mach_spec.P = 1;
    mach_spec.space_dimension = (MKL_INT)pow(2, mach_spec.num_qubits);

    optimization_spec_t opt_spec;
    opt_spec.ftol = 1e-16;
    opt_spec.xtol = 1e-16;
    opt_spec.nlopt_method = NLOPT_LN_NELDERMEAD;
    opt_spec.max_evals = 200 * mach_spec.P;

    cost_data_t cost_data;
    cost_data.cx_range = mach_spec.space_dimension;
    cost_data.x_range = mach_spec.space_dimension;
    cost_data.num_vertices = mach_spec.num_qubits;
    cost_data.graph = mkl_malloc(sizeof(MKL_INT) * mach_spec.num_qubits * mach_spec.num_qubits, DEF_ALIGNMENT);
    generate_graph(cost_data.graph, mach_spec.num_qubits, 0.5);
    print_graph(&cost_data, stdout);
    qaoa(&mach_spec, &cost_data, &opt_spec, &run_spec);
    mkl_free(cost_data.graph);
    return 0;
}
