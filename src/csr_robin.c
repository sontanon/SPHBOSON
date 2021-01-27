#include "tools.h"

// Finite difference coefficients for fourth order.
// One-sided.
#define O10 (+0.25)
#define O11 (-4.0 / 3.0)
#define O12 (+3.0)
#define O13 (-4.0)
#define O14 (25.0 / 12.0)

// Decay along r direction.
void robin_4th_order(
	double *aa,
	MKL_INT *ia,
	MKL_INT *ja,
	const MKL_INT offset,
	const MKL_INT dim,
	const MKL_INT ghost,
	const MKL_INT g_num,
	const MKL_INT i,
	const double dr,
	const MKL_INT n
)
{
	// Grid offset.
	MKL_INT k = g_num * dim;

	// Normalized coordinate values.
	double ri;

	// Row starts at offset.
	ia[k + i] = BASE + offset;

	// Coordinates.
	ri = (double)i + 0.5 - ghost;

	// Set values.
	aa[offset + 0] = O10 * ri;
	aa[offset + 1] = O11 * ri;
	aa[offset + 2] = O12 * ri;
	aa[offset + 3] = O13 * ri;
	aa[offset + 4] = O14 * ri + (double)n;

	// Column indices.
	ja[offset + 0] = BASE + k + i - 4;
	ja[offset + 1] = BASE + k + i - 3;
	ja[offset + 2] = BASE + k + i - 2;
	ja[offset + 3] = BASE + k + i - 1;
	ja[offset + 4] = BASE + k + i    ;

	// All done.
	return;
} 