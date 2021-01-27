#include "tools.h"

void symmetry(double *aa,
	MKL_INT *ia, 
	MKL_INT *ja, 
	const MKL_INT offset,
	const MKL_INT dim, 
	const MKL_INT g_num, 
	const MKL_INT i, 
	const MKL_INT ghost,
	const MKL_INT r_sym)
{
	// Row starts at offset.
	ia[g_num * dim + i] = BASE + offset;

	// Set values.
	aa[offset + 0] = 1.0;
	aa[offset + 1] = -(double)r_sym;

	// Set column values.
	ja[offset + 0] = BASE + g_num * dim + i;
	ja[offset + 1] = BASE + g_num * dim + 2 * ghost - (i + 1);

	// All done.
	return;
}
