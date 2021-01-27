#include "tools.h"
#include "omega_calc.h"

// Finite difference coefficients for fourth order.
// One-sided.
#define O10 (+0.25)
#define O11 (-4.0 / 3.0)
#define O12 (+3.0)
#define O13 (-4.0)
#define O14 (25.0 / 12.0)

// Decay along r direction.
void exp_decay_4th_order(
	double *aa,
	MKL_INT *ia,
	MKL_INT *ja,
	const MKL_INT offset,
	const MKL_INT dim,
	const MKL_INT ghost,
	const MKL_INT g_num,
	const MKL_INT i,
	const double dr,
	double *u,
	const MKL_INT w_idx,
	const double m
)
{
	// Grid offset.
	MKL_INT k = g_num * dim;

	// Normalized coordinate values.
	double ri, r;

	// Row starts at offset.
	ia[k + i] = BASE + offset;

	// Coordinates.
	ri = (double)i + 0.5 - ghost;
	r  = dr * ri;

	// Omega.
	double v = u[w_idx];
	double w = omega_calc(v, m);
	double w2 = w * w;
	double m2 = m * m;
	double chi = sqrt(m2 - w2);

	// RHS.
	double phi = u[k + i];
	double Dr_phi = O10 * u[k + i - 4] + O11 * u[k + i - 3] + O12 * u[k + i - 2] + O13 * u[k + i - 1] + O14 * phi;
	double f = ri * Dr_phi + (r * chi + 1.0) * phi;

	// Set values.
	aa[offset + 0] = O10 * ri;
	aa[offset + 1] = O11 * ri;
	aa[offset + 2] = O12 * ri;
	aa[offset + 3] = O13 * ri;
	aa[offset + 4] = O14 * ri + (1.0 + r * chi);
	aa[offset + 5] = dw_du(v, m) * (r * f + r * phi) * (-w / chi);

	// Column indices.
	ja[offset + 0] = BASE + k + i - 4;
	ja[offset + 1] = BASE + k + i - 3;
	ja[offset + 2] = BASE + k + i - 2;
	ja[offset + 3] = BASE + k + i - 1;
	ja[offset + 4] = BASE + k + i    ;
	ja[offset + 5] = BASE + w_idx;

	// All done.
	return;
} 