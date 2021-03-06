// ARCHITECTURE
#include "arch.h"

// Standard headers.
#include <stdio.h>
#include <stdlib.h>
#ifdef WIN
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <string.h>
#include <errno.h>

// System headers.
#include <sys/types.h>
#include <sys/stat.h>
#ifdef WIN
#include <io.h>
#include <direct.h>
#else
#include <unistd.h>
#endif

// Intel MKL
#include "mkl.h"

// Libconfig for parameter parsing.
#include <libconfig.h>

// MIN/MAX macros.
#define MIN(X, Y) ((X) < (Y)) ? (X) : (Y)
#define MAX(X, Y) ((X) > (Y)) ? (X) : (Y)

// ABS macro.
#define ABS(X) ((X) < 0) ? -(X) : (X)

// CSR matrix index base.
#define BASE 1

/* Macro for array sum z = alpha * x + beta * y: for alpha, beta scalars; z, x, y arrays. */
#define ARRAY_SUM(Z, ALPHA, X, BETA, Y) array_sum((Z), (ALPHA), (X), (BETA), (Y), dim)
void array_sum(double *z, const double alpha, const double *x, const double beta, const double *y, const MKL_INT dim);

// Safe allocation macros.
#define SAFE_MALLOC(n) safe_malloc((n), __FILE__, __LINE__)
void *safe_malloc(const size_t n, const char *file, const MKL_INT line);

// Safe deallocation macros.
#define SAFE_FREE(x) safe_free((x), __FILE__, __LINE__)
void safe_free(void *x, const char *file, const MKL_INT line);

// CSR matrix type.
typedef struct csr_matrices
{
	double *a;
	MKL_INT *ia;
	MKL_INT *ja;
	MKL_INT nrows;
	MKL_INT ncols;
	MKL_INT nnz;
	MKL_INT analysis_phase;

} csr_matrix;

// Forward declarations.
// 
// Write simple ASCII 1D file.
void write_single_file(const double *u, const char *fname, const MKL_INT dim);
void write_single_file_integer(const MKL_INT *u, const char *fname, const MKL_INT dim);
// Read simple ASCII 1D file.
void read_single_file(double *u, const char *fname, const MKL_INT dim, const char *source_file, const MKL_INT source_line);
// Create CSR matrix.
void csr_allocate(csr_matrix *A, const MKL_INT nrows, const MKL_INT ncols, const MKL_INT nnz);
// Destroy CSR matrx.
void csr_deallocate(csr_matrix *A);
// CSR matrix print.
void csr_print(csr_matrix *A, const char *vA, const char *iA, const char *jA);
