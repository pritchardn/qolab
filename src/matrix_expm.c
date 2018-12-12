//
// Created by nicholas on 6/07/18.
//

#include "matrix_expm.h"
#include "globals.h"

void spmatrix_expm_z_diag(const MKL_Complex16 *diag, double alpha, MKL_INT nnz, MKL_Complex16 *state) {
    MKL_Complex16 *tempValues = mkl_calloc((size_t) nnz, sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(tempValues);
    MKL_Complex16 *resultValues = mkl_malloc(nnz * sizeof(MKL_Complex16), DEF_ALIGNMENT);
    check_alloc(resultValues);

    cblas_zcopy(nnz, diag, 1, tempValues, 1);
    cblas_zdscal(nnz, alpha, tempValues, 1);
    vzExp(nnz, tempValues, resultValues);
    cblas_zcopy(nnz, state, 1, tempValues, 1);
    vzMul(nnz, resultValues, tempValues, state);

    mkl_free(tempValues);
    mkl_free(resultValues);
}

void spmatrix_expm_cheby(sparse_matrix_t *matrix, MKL_Complex16 *state, MKL_Complex16 dt,
                         MKL_Complex16 minE, MKL_Complex16 maxE,
                         MKL_INT side_len) {
    int i, terms;
    double alpha;
    complex double emin, emax, t, EmEm, d2EmEm, imagM, neg1, bessj0, bessj1, bessjn, ztemp;
    sparse_status_t status;
    MKL_Complex16 mkl_ztemp1, mkl_ztemp2;
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    MKL_Complex16 **work = mkl_malloc(4 * sizeof(MKL_Complex16 *), DEF_ALIGNMENT);
    check_alloc(work);

    for (i = 0; i < 4; ++i) {
        work[i] = mkl_calloc((size_t) side_len, sizeof(MKL_Complex16), DEF_ALIGNMENT);
        check_alloc(work[i]);
    }

    emin = minE.real + I * minE.imag;
    emax = maxE.real + I * maxE.imag;
    t = dt.real + I * dt.imag;

    EmEm = (emax + emin) / (emax - emin);
    d2EmEm = -2.0 / (emax - emin);
    alpha = creal(I * (emax - emin) * t / 2.0);

    cblas_zcopy(side_len, state, 1, work[0], 1);

    status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, (MKL_Complex16) {1.0, 0.0}, *matrix, descr, work[0],
                             (MKL_Complex16) {0.0, 0.0}, work[1]);
    mkl_error_parse(status, stderr);

    mkl_ztemp1.real = creal(EmEm);
    mkl_ztemp1.imag = cimag(EmEm);

    mkl_ztemp2.real = creal(d2EmEm);
    mkl_ztemp2.imag = cimag(d2EmEm);

    cblas_zaxpby(side_len, &mkl_ztemp1, work[0], 1, &mkl_ztemp2, work[1], 1);

    bessj0 = jn(0, alpha);
    bessj1 = jn(1, alpha);
    bessj1 *= 2 * I;

    mkl_ztemp1.real = creal(bessj1);
    mkl_ztemp1.imag = cimag(bessj1);

    mkl_ztemp2.real = creal(bessj0);
    mkl_ztemp2.imag = cimag(bessj0);

    cblas_zcopy(side_len, state, 1, work[3], 1);

    cblas_zaxpby(side_len, &mkl_ztemp1, work[1], 1, &mkl_ztemp2, work[3], 1);

    terms = 0;
    while (fabs(2.0 * jn(terms, alpha)) > 1e-18) {
        terms++;
    }

    EmEm *= 2.0;
    d2EmEm *= 2.0;
    imagM = 2 * I * I;
    neg1 = -1.0;


    for (i = 2; i <= terms; ++i) {
        status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, (MKL_Complex16) {1.0, 0.0}, *matrix, descr, work[1],
                                 (MKL_Complex16) {0.0, 0.0}, work[2]);
        mkl_error_parse(status, stderr);

        mkl_ztemp1.real = creal(EmEm);
        mkl_ztemp1.imag = cimag(EmEm);

        mkl_ztemp2.real = creal(d2EmEm);
        mkl_ztemp2.imag = cimag(d2EmEm);

        cblas_zaxpby(side_len, &mkl_ztemp1, work[1], 1, &mkl_ztemp2, work[2], 1);

        mkl_ztemp1.real = creal(neg1);
        mkl_ztemp1.imag = cimag(neg1);

        cblas_zaxpy(side_len, &mkl_ztemp1, work[0], 1, work[2], 1);

        bessjn = jn(i, alpha);

        ztemp = imagM * bessjn;
        mkl_ztemp1.real = creal(ztemp);
        mkl_ztemp1.imag = cimag(ztemp);

        cblas_zaxpy(side_len, &mkl_ztemp1, work[2], 1, work[3], 1);

        imagM *= I;

        cblas_zcopy(side_len, work[1], 1, work[0], 1);
        cblas_zcopy(side_len, work[2], 1, work[1], 1);
    }

    ztemp = cexp(-I * (emax + emin) * (t / 2.0));

    mkl_ztemp1.real = creal(ztemp);
    mkl_ztemp1.imag = cimag(ztemp);

    cblas_zscal(side_len, &mkl_ztemp1, work[3], 1);
    cblas_zcopy(side_len, work[3], 1, state, 1);


    mkl_free(work[3]);
    mkl_free(work[2]);
    mkl_free(work[1]);
    mkl_free(work[0]);
    mkl_free(work);
}
