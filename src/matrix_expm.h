//
// Created by nicholas on 6/07/18.
//

#ifndef GRAPHSIMILARITY_MATRIX_EXPM_H
#define GRAPHSIMILARITY_MATRIX_EXPM_H

#include <mkl.h>
#include <mathimf.h>
#include <complex.h>

void spmatrix_expm_z_diag(const MKL_Complex16 *diag, double alpha, MKL_INT nnz, MKL_Complex16 *state);

void spmatrix_expm_cheby(sparse_matrix_t *matrix, MKL_Complex16 *state, MKL_Complex16 dt,
                         MKL_Complex16 minE, MKL_Complex16 maxE, int side_len);

#endif //GRAPHSIMILARITY_MATRIX_EXPM_H
