/**
 * @author Nicholas Pritchard
 * @date 22/03/2019
 */

#ifndef QOLAB_EIGEN_SOLVE_H
#define QOLAB_EIGEN_SOLVE_H

#include <mkl.h>

void test_eigen_solve(sparse_matrix_t *A, MKL_INT rows, MKL_INT cols, MKL_INT nnz);
#endif //QOLAB_EIGEN_SOLVE_H
