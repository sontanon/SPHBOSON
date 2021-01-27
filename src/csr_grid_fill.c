#include "tools.h"
#include "csr_symmetry.h"
#include "csr_robin.h"
#include "csr_exp_decay.h"

#define GNUM 3

#define EVEN	 1
#define ODD	-1

#define EXP_DECAY_TYPE 	0
#define ROBIN_TYPE_1 	1

void csr_grid_fill_4th(
	csr_matrix A, 
	const MKL_INT dim,
	const MKL_INT ghost,
	const double dr,
	double *u,
	const double m,
	const MKL_INT r_sym[GNUM],
	const MKL_INT bound_order[GNUM],
	const MKL_INT nnz[GNUM],
	const MKL_INT p_c[GNUM],
	const MKL_INT p_s[GNUM],
	const MKL_INT p_bound[GNUM],
	void (*j_c)(double *, MKL_INT *, MKL_INT *,
		const MKL_INT, const MKL_INT, const MKL_INT, const double,
		const double, const double,
		double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double,
		const MKL_INT, const MKL_INT, const MKL_INT),
	void (*j_s)(double *, MKL_INT *, MKL_INT *,
		const MKL_INT, const MKL_INT, const MKL_INT, const double,
		const double, const double,
		double, double, double, double, double, double,
		double, double, double, double, double, double,
		double, double, double, double, double, double,
		const MKL_INT, const MKL_INT, const MKL_INT)
)
{
	// Auxiliary integers.
	MKL_INT i = 0, k = 0;

	// Omega index.
	MKL_INT w_idx = GNUM * dim;

	// Integer arrays for offsets and start indices.
	MKL_INT offset[GNUM] = { 0, 0, 0 };
	MKL_INT t_offset[GNUM] = { 0, 0, 0 };
	MKL_INT start_offset[GNUM] = { 0, 0, 0 };

	// Set start/initial offsets.
	offset[0] = start_offset[0] = 0;
	for (k = 1; k < GNUM; ++k)
	{
		offset[k] = start_offset[k] = offset[k - 1] + nnz[k - 1];
	}

	// Fetch xi variable once.
	double xi = u[w_idx];

	// Axis symmetry.
	for (i = 0; i < ghost; ++i)
	{
		for (k = 0; k < GNUM; ++k)
		{
			// Increase offsets by 2.
			symmetry(A.a, A.ia, A.ja, offset[k], dim, k, i, ghost, r_sym[k]);
			offset[k] += 2;
		}
	}
	// Set temporary offset.
	for (k = 0; k < GNUM; ++k)
		t_offset[k] = offset[k];

	// Fill interior points.
	#pragma omp parallel shared(A) private(offset, i, k)
	{
		#pragma omp for schedule(dynamic, 1)
		for (i = ghost; i < dim - 2; ++i)
		{
			// Each iteration of i loop will fill p_c values.
			for (k = 0; k < GNUM; ++k)
				offset[k] = t_offset[k] + (i - ghost) * p_c[k];

			// Fill matrix coefficients.
			(*j_c)(A.a, A.ia, A.ja,
				dim, ghost, i, dr, m, xi,
				u[          i - 2], u[          i - 1], u[          i], u[          i + 1], u[          i + 2],
				u[    dim + i - 2], u[    dim + i - 1], u[    dim + i], u[    dim + i + 1], u[    dim + i + 2],
				u[2 * dim + i - 2], u[2 * dim + i - 1], u[2 * dim + i], u[2 * dim + i + 1], u[2 * dim + i + 2],
				offset[0], offset[1], offset[2]);
		}
	}
	// At this point we have filled:
	for (k = 0; k < GNUM; ++k)
		offset[k] = start_offset[k] + 2 * ghost + (dim - 2 - ghost) * p_c[k];

	// Semi-onesided stencil for i = dim - 2.
	i = dim - 2;
	(*j_s)(A.a, A.ia, A.ja,
		dim, ghost, i, dr, m, xi,
		u[0 * dim + i - 4], u[0 * dim + i - 3], u[0 * dim + i - 2], u[0 * dim + i - 1], u[0 * dim + i], u[0 * dim + i + 1], 
		u[1 * dim + i - 4], u[1 * dim + i - 3], u[1 * dim + i - 2], u[1 * dim + i - 1], u[1 * dim + i], u[1 * dim + i + 1], 
		u[2 * dim + i - 4], u[2 * dim + i - 3], u[2 * dim + i - 2], u[2 * dim + i - 1], u[2 * dim + i], u[2 * dim + i + 1], 
		offset[0], offset[1], offset[2]);

	// Increase offset.
	for (k = 0; k < GNUM; ++k)
		offset[k] += p_s[k];

	// Boundary condition for i = dim - 1.
	i = dim - 1;
	robin_4th_order(A.a, A.ia, A.ja, offset[0], dim, ghost, 0, i, dr, bound_order[0]);
	robin_4th_order(A.a, A.ia, A.ja, offset[1], dim, ghost, 1, i, dr, bound_order[1]);
	exp_decay_4th_order(A.a, A.ia, A.ja, offset[2], dim, ghost, 2, i, dr, u, w_idx, m);

	// Increase offsets.
	for (k = 0; k < GNUM; ++k)
		offset[k] += p_bound[k];

	// All done.
	return;
}