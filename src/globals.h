#include <stdio.h>
#include <mkl.h>

#define DEF_ALIGNMENT 64


//Mathematical
MKL_INT factorial(int n);

//Administrative
void mkl_error_parse(int error, FILE *stream);
void check_alloc(void *pointer);

//Debugging
void mkl_sparse_print(sparse_matrix_t *matrix, FILE *stream);
