#include "globals.h"

MKL_INT factorial(int n) {
    MKL_INT result = 0;
    if (n == 0) {
        return 1;
    }
    if (n > 0) {
        result++;
        for (int i = 1; i <= n; ++i) {
            result *= i;
        }
    }
    return result;
}

void check_alloc(void *pointer) {
    if (pointer == NULL) {
        fprintf(stderr, "Alloc fail\n");
        exit(EXIT_FAILURE);
    }
}

void mkl_error_parse(int error_code, FILE *stream) {
    switch (error_code) {
        case SPARSE_STATUS_SUCCESS:
            //fprintf(stream, "%d SUCCESS\n", error_code);
            break;
        case SPARSE_STATUS_NOT_INITIALIZED:
            fprintf(stream, "%d NOT INITIALISED\n", error_code);
            exit(EXIT_FAILURE);
        case SPARSE_STATUS_ALLOC_FAILED:
            fprintf(stream, "%d ALLOC FAILED\n", error_code);
            exit(EXIT_FAILURE);
        case SPARSE_STATUS_INVALID_VALUE:
            fprintf(stream, "%d INVALID VALUE\n", error_code);
            exit(EXIT_FAILURE);
        case SPARSE_STATUS_EXECUTION_FAILED:
            fprintf(stream, "%d EXEC FAILED\n", error_code);
            exit(EXIT_FAILURE);
        case SPARSE_STATUS_INTERNAL_ERROR:
            fprintf(stream, "%d INTERNAL ERROR\n", error_code);
            exit(EXIT_FAILURE);
        case SPARSE_STATUS_NOT_SUPPORTED:
            fprintf(stream, "%d NOT SUPPORTED\n", error_code);
            exit(EXIT_FAILURE);
        default:
            fprintf(stream, "Something else\n");
            break;
    }
}