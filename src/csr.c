#include "tools.h"
#include "param.h"
#include "csr_grid_fill.h"
#include "csr_vars.h"
#include "csr_omega_constraint.h"

#define EVEN	 1
#define ODD	-1

#define ROBIN_TYPE_1 	1
#define EXP_DECAY_TYPE 	0

void nnz_jacobian_get_nnzs(MKL_INT *p_nnz1, MKL_INT *p_nnz2, MKL_INT *p_nnz3)
{
	MKL_INT nnz1, nnz2, nnz3;
	// Select order.
	if (order == 4)
	{
		// Interior points plus parity boundaries.
		nnz1 = 12 * NrInterior + 13 + 5 + 2 + 2;
		nnz2 = 12 * NrInterior + 13 + 5 + 2 + 2;
		nnz3 = 16 * NrInterior + 17 + 6 + 2 + 2;
	}
	*p_nnz1 = nnz1;
	*p_nnz2 = nnz2;
	*p_nnz3 = nnz3;

	return;
}

MKL_INT nnz_jacobian(void)
{
	// Number of nonzero elements per grid function.
	MKL_INT nnz1, nnz2, nnz3;

	nnz_jacobian_get_nnzs(&nnz1, &nnz2, &nnz3);

	// Total number of nonzeros includes one extra for omega.
	return nnz1 + nnz2 + nnz3 + 1;
}

void csr_gen_jacobian(csr_matrix A, double *u, const MKL_INT print)
{	
	// Number of elements we have filled in.
	MKL_INT offset = 0;

	// Number of nonzero elements per grid functions.
	MKL_INT nnz1, nnz2, nnz3;

	// Calculate nonzeros.
	nnz_jacobian_get_nnzs(&nnz1, &nnz2, &nnz3);

	// Integer arrays.
	MKL_INT r_sym[GNUM] = { EVEN, EVEN, EVEN };
	MKL_INT bound_order[GNUM] = { ROBIN_TYPE_1, ROBIN_TYPE_1, EXP_DECAY_TYPE };
	MKL_INT nnzs[GNUM] = { nnz1, nnz2, nnz3 };

	MKL_INT p_c[GNUM] = { 12, 12, 16 };
	MKL_INT p_s[GNUM] = { 13, 13, 17 };
	MKL_INT p_bound[GNUM] = { 5, 5, 6 };

	if (order == 4)
	{
		csr_grid_fill_4th(A, dim, ghost, dr, u, m, r_sym, bound_order, nnzs, p_c, p_s, p_bound, jacobian_4th_order_variable_omega_c, jacobian_4th_order_variable_omega_s);
	}

	// FINALLY, FILL OMEGA EQUATION OR u3(fixed_phi) CONSTRAINT.
	offset = nnz1 + nnz2 + nnz3;
	omega_constraint(A.a, A.ia, A.ja, offset, dim, GNUM, w_idx, fixedPhi, ghost);

	// FILL LAST ELEMENT WITH TOTAL NUMBER OF NONZEROS.
	A.ia[w_idx + 1] = BASE + A.nnz;

	if (print)
		csr_print(&A, "a.asc", "ia.asc", "ja.asc");

	// All done.
	return;
}
