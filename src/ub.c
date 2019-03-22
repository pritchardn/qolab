/**
 * @author Nicholas Pritchard
 * @date 1/07/2018
 */
#include "ub.h"

/**
 * @brief Generates the driver hamiltonian for a given problem with double values.
 * @details Defines the continuous time quantum walk that allows for 'probability' to flow around candidate solution bitstrings
 * In the standard QAOA this defines a fully connected hyper-cube, in a restricted QAOA this is a problem dependent
 * subset of this graph. The double values allow us to quickly solve for the eigenvalues of this matrix
 * @param meta_data Describes the full simulation. num_qubits, cost_data are used
 * @param mask (optional) Returns true given a valid input, false otherwise.
 * @warning Will need to be converted to a complex matrix before use with the QAOA module
 */
MKL_INT generate_ub(qaoa_data_t *meta_data, bool (*mask)(unsigned int, cost_data_t *cost_data)) {
    sparse_status_t status;
    MKL_INT nnz = 0;
    MKL_INT space_dimension = meta_data->machine_spec->space_dimension;
    double *values = mkl_calloc((size_t) space_dimension * (meta_data->machine_spec->num_qubits + 1),
                                sizeof(double),
                                DEF_ALIGNMENT);
    MKL_INT *row_begin = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *row_end = mkl_calloc((size_t) space_dimension + 1, sizeof(MKL_INT), DEF_ALIGNMENT);
    MKL_INT *col_index = mkl_calloc((size_t) space_dimension * (meta_data->machine_spec->num_qubits) + 1,
                                    sizeof(MKL_INT), DEF_ALIGNMENT);
    check_alloc(values);
    check_alloc(row_begin);
    check_alloc(row_end);
    check_alloc(col_index);
    for (unsigned int i = 0; i < space_dimension; ++i) {
        row_begin[i] = nnz;
        for (unsigned int j = 0; j < meta_data->machine_spec->num_qubits; ++j) {
            unsigned int col = ((i ^ ((unsigned int) 1 << j)));
            if (mask(col, meta_data->cost_data)) {
                values[nnz] = 1.0;
                row_end[i] = nnz;
                col_index[nnz] = (MKL_INT) col;
                nnz++;
            }
        }
    }
    status = mkl_sparse_d_create_csr(&meta_data->ub, (sparse_index_base_t) SPARSE_INDEX_BASE_ZERO, \
    space_dimension, space_dimension, row_begin, row_end, col_index, values);
    mkl_error_parse(status, stdout);

    return nnz;
}

/**
 * @brief Converts a double valued UB matrix to a complex valued UB matrix
 * @param meta_data The structure which holds the UB matrix itself
 * @param ub_nnz The number of non-zero elements in the UB matrix
 */
void convert_ub(qaoa_data_t *meta_data, MKL_INT ub_nnz) {
    sparse_status_t status;
    sparse_index_base_t index_base;
    MKL_INT rows;
    MKL_INT cols;
    MKL_INT *rows_start;
    MKL_INT *rows_end;
    MKL_INT *col_indx;
    double *values;
    MKL_Complex16 *new_values = mkl_calloc((size_t) ub_nnz, sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(new_values);

    status = mkl_sparse_d_export_csr(meta_data->ub, &index_base, &rows, &cols, &rows_start, &rows_end, &col_indx,
                                     &values);
    mkl_error_parse(status, stderr);

    for (MKL_INT i = 0; i < ub_nnz; ++i) {
        //Shortcut to avoid * -I later
        new_values[i].imag = -1.0;
    }

    status = mkl_sparse_z_create_csr(&meta_data->ub, index_base, rows, cols, rows_start, rows_end, col_indx,
                                     new_values);
    mkl_error_parse(status, stderr);

    mkl_free(values);
}