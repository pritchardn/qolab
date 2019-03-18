#include "globals.h"

/**
 * Helper function which returns the factorial of an int.
 * TODO: Unit test for overflow
 * @param n The integer have the factorial computed on
 * @return MKL_INT The factorial of n
 * @warning This is prone to overflow
 */
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

/**
 * Helper function to check the validity of a pointer and failes gracefully
 * @param pointer The memory to be checked
 */
void check_alloc(void *pointer) {
    if (pointer == NULL) {
        fprintf(stderr, "Allocation fail\n");
        exit(EXIT_FAILURE);
    }
}

/**
 * A custom mkl error code parser.
 *
 * Checks for a variety of possible errors and exits gracefully:
 *   - SUCCESS
 *   - NOT INITIALISED
 *   - ALLOC FAILED
 *   - INVALID VALUE
 *   - EXEC_FAILED
 *   - INTERNAL ERROR
 *   - NOT SUPPORTED
 * @param error_code The MKL error flag to be processed
 * @param stream The specified stream to send debug messages to
 * TODO: Test for file, stream and NULL
 */
void mkl_error_parse(int error_code, FILE *stream) {
    if (stream == NULL) {
        stream = stderr;
    }
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