/**
 * @author Nicholas Pritchard
 * @date 22/03/2019
 */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

void max_eigen_find() {
/* Matrix A of size N in CSR format */
    MKL_INT N = 4;               /* number of rows in matrix A */
    MKL_INT M = 4;               /* number of columns in matrix A */
    MKL_INT nnz = 8;             /* number of non-zeros in matrix */

    MKL_INT ia[5] = {1, 3, 5, 7, 9};                         /* ia array from CSR format */
    MKL_INT ja[8] = {1, 2, 1, 2, 3, 4, 3, 4};                   /* ja array from CSR format */
    double a[8] = {6.0, 2.0, 2.0, 3.0, 2.0, -1.0, -1.0, 2.0}; /* val array from CSR format */

    double Eig[4] = {1.0, 2.0, 3.0, 7.0}; /* Exact eigenvalues */

    /* mkl_sparse_d_ev input parameters */
    char which = 'S'; /* Which eigenvalues to calculate. ('L' - largest (algebraic) eigenvalues, 'S' - smallest (algebraic) eigenvalues) */
    MKL_INT pm[128];     /* This array is used to pass various parameters to Extended Eigensolver Extensions routines. */
    MKL_INT k0 = 3;     /* Desired number of max/min eigenvalues */

    /* mkl_sparse_d_ev output parameters */
    MKL_INT k;           /* Number of eigenvalues found (might be less than k0). */
    double E[4];   /* Eigenvalues */
    double X[4];   /* Eigenvectors */
    double res[4]; /* Residual */

    /* Local variables */
    MKL_INT info;               /* Errors */
    MKL_INT compute_vectors = 0;/* Flag to compute eigenvecors */
    MKL_INT tol = 7;            /* Tolerance */
    double Y[4];               /* Y=(X')*X-I */
    double sparsity;           /* Sparsity of randomly generated matrix */
    int i, j;
    double smax, t;

    /* Input variables for DGEMM */
    char DGEMMC = 'T';       /* Character for GEMM routine, transposed case */
    char DGEMMN = 'N';       /* Character for GEMM routine, non-transposed case */
    double one = 1.0;         /* alpha parameter for GEMM */
    double zero = 0.0;         /* beta  parameter for GEMM */
    MKL_INT ldx = N;           /* Leading dimension for source arrays in GEMM */
    MKL_INT ldy;                /* Leading dimension for destination array in GEMM */

    /* Sparse BLAS IE variables */
    sparse_status_t status;
    sparse_matrix_t A = NULL; /* Handle containing sparse matrix in internal data structure */
    struct matrix_descr descr; /* Structure specifying sparse matrix properties */

    /* Create handle for matrix A stored in CSR format */
    descr.type = SPARSE_MATRIX_TYPE_GENERAL; /* Full matrix is stored */
    status = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ONE, N, N, ia, ia + 1, ja, a);
    if (status != 0) {
        printf("Building failed\n");
        exit(EXIT_FAILURE);
    }
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

    mkl_sparse_destroy(A);
}
