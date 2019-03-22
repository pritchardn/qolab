/**
 * @author Nicholas Pritchard
 * @date 22/03/2019
 */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "globals.h"

double max_eigen_find(sparse_matrix_t A) {
    double result;
    double Eig[4] = {1.0, 2.0, 3.0, 7.0}; /* Exact eigenvalues */

    MKL_INT rows;
    MKL_INT cols;
    MKL_INT *rows_start;
    MKL_INT *rows_end;
    MKL_INT *col_indx;
    double *values;
    sparse_status_t status;
    sparse_index_base_t index_base;

    status = mkl_sparse_d_export_csr(A, &index_base, &rows, &cols, &rows_start, &rows_end, &col_indx,
                                     &values);
    mkl_error_parse(status, stderr);

    /* mkl_sparse_d_ev input parameters */
    char which = 'S'; /* Which eigenvalues to calculate. ('L' - largest (algebraic) eigenvalues, 'S' - smallest (algebraic) eigenvalues) */
    MKL_INT pm[128];     /* This array is used to pass various parameters to Extended Eigensolver Extensions routines. */
    MKL_INT k0 = 3;     /* Desired number of max/min eigenvalues */

    /* mkl_sparse_d_ev output parameters */
    MKL_INT k;           /* Number of eigenvalues found (might be less than k0). */
    double *E = mkl_calloc((size_t) rows, sizeof(double), DEF_ALIGNMENT);   /* Eigenvalues */
    double *X = mkl_calloc((size_t) rows, sizeof(double), DEF_ALIGNMENT);   /* Eigenvectors */
    double *res = mkl_calloc((size_t) rows, sizeof(double), DEF_ALIGNMENT); /* Residual */
    check_alloc(E);
    check_alloc(X);
    check_alloc(res);

    /* Local variables */
    MKL_INT info;               /* Errors */
    MKL_INT compute_vectors = 0;/* Flag to compute eigenvecors */
    MKL_INT tol = 7;            /* Tolerance */
    int i;

    /* Sparse BLAS IE variables */
    struct matrix_descr descr; /* Structure specifying sparse matrix properties */

    /* Create handle for matrix A stored in CSR format */
    descr.type = SPARSE_MATRIX_TYPE_GENERAL; /* Full matrix is stored */


    /* Step 2. Call mkl_sparse_ee_init to define default input values */
    mkl_sparse_ee_init(pm);

    pm[1] = tol; /* Set tolerance */
    pm[6] = compute_vectors;

    /* Step 3. Solve the standard Ax = ex eigenvalue problem. */
    info = mkl_sparse_d_ev(&which, pm, A, descr, k0, &k, E, X, res);

    printf("mkl_sparse_d_ev output info %d \n", info);
    if (info != 0) {
        printf("Routine mkl_sparse_d_ev returns code of ERROR: %i", (int) info);
        exit(EXIT_FAILURE);
    }
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("#mode found/subspace %d %d \n", k, k0);
    printf("Index/Exact Eigenvalues/Estimated Eigenvalues/Residuals\n");
    for (i = 0; i < k; i++) {
        printf("   %d  %.15e %.15e %.15e \n", i, Eig[i], E[i], res[i]);
    }

    result = E[0];

    mkl_free(E);
    mkl_free(X);
    mkl_free(res);

    return result;
}
