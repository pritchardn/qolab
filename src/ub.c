//
// Created by nicholas on 1/07/18.
//

#include "ub.h"

/**
 * Generates initial state defining a quantum walk on our qubits. Assumes all qubits are fully connected, additional
 * walk masks can be applied (coming later). Builds matrix in coordinate format and then converts to CSC format
 * @param graph_size
 * @param ub
 */
void generate_ub(int num_qubits, sparse_matrix_t *ub) {
    sparse_status_t status;
    MKL_INT nnz = 0;
    MKL_INT space_dimension = (MKL_INT) pow(2, num_qubits);
    MKL_Complex16 *values = mkl_calloc((size_t) space_dimension * (num_qubits + 1), sizeof(MKL_Complex16),
                                       DEF_ALIGNMENT);
    MKL_INT *col_begin = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *col_end = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *row_index = mkl_calloc((size_t) space_dimension * (num_qubits) + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(values);
    check_alloc(col_begin);
    check_alloc(col_end);
    check_alloc(row_index);
    for (unsigned int i = 0; i < space_dimension; ++i) {
        col_begin[i] = nnz;
        for (unsigned int j = 0; j < num_qubits; ++j) {
            unsigned int row = ((i ^ ((unsigned int) 1 << j)));
            values[nnz].imag = -1.0;
            col_end[i] = nnz;
            row_index[nnz] = (MKL_INT) row;
            nnz++;
        }
    }
    status = mkl_sparse_z_create_csc(ub, (sparse_index_base_t) SPARSE_INDEX_BASE_ZERO, \
    space_dimension, space_dimension, col_begin, col_end, row_index, values);
    mkl_error_parse(status, stdout);
}


void generate_ub_mpi(int graph_size, sparse_matrix_t *ub) {
    sparse_status_t status;
    MKL_INT nnz = 0;
    MKL_INT num_qubits = (MKL_INT) ceill(log2l(factorial(graph_size)));
    MKL_INT space_dimension = (MKL_INT) pow(2, num_qubits);
    MKL_Complex16 *values = mkl_calloc((size_t) space_dimension * (num_qubits + 1), sizeof(MKL_Complex16),
                                       DEF_ALIGNMENT);
    MKL_INT *col_begin = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *col_end = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *row_index = mkl_calloc((size_t) space_dimension * (num_qubits) + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(values);
    check_alloc(col_begin);
    check_alloc(col_end);
    check_alloc(row_index);
    for (unsigned int i = 0; i < space_dimension; ++i) {
        col_begin[i] = nnz;
        for (unsigned int j = 0; j < num_qubits; ++j) {
            unsigned int row = ((i ^ ((unsigned int) 1 << j)));
            values[nnz].imag = -1.0;
            col_end[i] = nnz;
            row_index[nnz] = (MKL_INT) row;
            nnz++;
        }
    }
    status = mkl_sparse_z_create_csc(ub, (sparse_index_base_t) SPARSE_INDEX_BASE_ZERO, \
    space_dimension, space_dimension, col_begin, col_end, row_index, values);
    mkl_error_parse(status, stdout);
}
