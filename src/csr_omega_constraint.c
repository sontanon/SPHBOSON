#include "tools.h"

// Fill-in omega constraint.
void omega_constraint(double *aa, 
	MKL_INT *ia, 
	MKL_INT *ja,
	const MKL_INT offset, 
	const MKL_INT dim,
	const MKL_INT g_num,
	const MKL_INT w_idx,
	const MKL_INT fixedPhi,
	const MKL_INT i)
{
	ia[w_idx]  = BASE + offset;
	aa[offset] = 1.0;
	ja[offset] = (fixedPhi) ? BASE + (g_num - 1) * dim + i : BASE + w_idx;

	return;
}
