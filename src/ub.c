/**
 * @author Nicholas Pritchard
 * @date 1/07/2018
 */
#include "ub.h"

/**
 * @brief Generates the driver hamiltonian for a given problem.
 * @details Defines the continuous time quantum walk that allows for 'probability' to flow around candidate solution bitstrings
 * In the standard QAOA this defines a fully connected hyper-cube, in a restricted QAOA this is a problem dependent
 * subset of this graph
 * @param meta_data Describes the full simulation. num_qubits, cost_data are used
 * @param mask (optional) Returns true given a valid input, false otherwise.
 */
void generate_ub(qaoa_data_t *meta_data, bool (*mask)(unsigned int, cost_data_t *cost_data)) {
    sparse_status_t status;
    MKL_INT nnz = 0;
    MKL_INT space_dimension = meta_data->machine_spec->space_dimension;
    MKL_Complex16 *values = mkl_calloc((size_t) space_dimension * (meta_data->machine_spec->num_qubits + 1),
                                       sizeof(MKL_Complex16),
                                       DEF_ALIGNMENT);
    MKL_INT *col_begin = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *col_end = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *row_index = mkl_calloc((size_t) space_dimension * (meta_data->machine_spec->num_qubits) + 1,
                                    sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(values);
    check_alloc(col_begin);
    check_alloc(col_end);
    check_alloc(row_index);
    for (unsigned int i = 0; i < space_dimension; ++i) {
        col_begin[i] = nnz;
        for (unsigned int j = 0; j < meta_data->machine_spec->num_qubits; ++j) {
            unsigned int row = ((i ^ ((unsigned int) 1 << j)));
            if (mask(row, meta_data->cost_data)) {
                values[nnz].imag = -1.0;
                col_end[i] = nnz;
                row_index[nnz] = (MKL_INT) row;
                nnz++;
            }
        }
    }
    status = mkl_sparse_z_create_csc(&meta_data->ub, (sparse_index_base_t) SPARSE_INDEX_BASE_ZERO, \
    space_dimension, space_dimension, col_begin, col_end, row_index, values);
    mkl_error_parse(status, stdout);
}