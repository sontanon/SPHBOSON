#include "tools.h"
#include "param.h"

#include "derivatives.h"
#include "omega_calc.h"

void rhs(double *f, double *u)
{
	// Omega.
	double w = omega_calc(u[w_idx], m);
	double w2 = w * w;
	double m2 = m * m;

	// Auxiliary doubles.
	double r;
	double l_alpha, Dr_l_alpha, Drr_l_alpha, alpha2;
	double l_psi, Dr_l_psi, Drr_l_psi, psi4;
	double phi, Dr_phi, Drr_phi, phi2;
	//double scale;

	// Loop counters.
	MKL_INT i = 0, k = 0;

	// Calculate derivatives.
	diff1r(Dr_u          , u          , 1);
	diff1r(Dr_u +     dim, u +     dim, 1);
	diff1r(Dr_u + 2 * dim, u + 2 * dim, 1);

	diff2r(Drr_u          , u          , 1);
	diff2r(Drr_u +     dim, u +     dim, 1);
	diff2r(Drr_u + 2 * dim, u + 2 * dim, 1);

	// Parity on axis.
	for (i = 0; i < ghost; ++i)
	{
		for (k = 0; k < GNUM; ++k)
			f[k * dim + i] = 0.0;
	}

	// Main interior points.
	#pragma omp parallel shared(f) private(r,\
		    l_alpha, Dr_l_alpha, Drr_l_alpha,\
		    l_psi, Dr_l_psi, Drr_l_psi,\
		    phi, Dr_phi, Drr_phi,\
		    alpha2, psi4, phi2)
	{
		#pragma omp for schedule(dynamic, 1)
		for (i = ghost; i < dim - 1; ++i)
		{
			// Get coordinate value.
			r = dr * ((double)i + 0.5 - ghost);

			// Fetch values.
			l_alpha     =     u[i];
			Dr_l_alpha  =  Dr_u[i];
			Drr_l_alpha = Drr_u[i];
			l_psi       =     u[dim + i];
			Dr_l_psi    =  Dr_u[dim + i];
			Drr_l_psi   = Drr_u[dim + i];
			phi         =     u[2 * dim + i];
			Dr_phi      =  Dr_u[2 * dim + i];
			Drr_phi     = Drr_u[2 * dim + i];
			// Common values.
			alpha2 = exp(2.0 * l_alpha);
			psi4   = exp(4.0 * l_psi);
			phi2   = phi * phi;

			// Calculate RHS's.
			f[          i] = -dr * dr * (Drr_l_alpha + 2.0 * (Dr_l_alpha / r) 
					+ Dr_l_alpha * (Dr_l_alpha + 2.0 * Dr_l_psi) 
					+ 4.0 * M_PI * psi4 * (m2 - 2.0 * w2 / alpha2) * phi2);
			f[    dim + i] = -dr * dr * (Drr_l_psi + 2.0 * (Dr_l_psi / r) 
					+ Dr_l_psi * Dr_l_psi 
					+ M_PI * (Dr_phi * Dr_phi + psi4 * (m2 + w2 / alpha2) * phi2));
			f[2 * dim + i] = -dr * dr * (Drr_phi + 2.0 * (Dr_phi / r) 
					+ Dr_phi * (Dr_l_alpha + 2.0 * Dr_l_psi) 
					+ psi4 * (-m2 + w2 / alpha2) * phi);

		}
	}

	// Robin boundary condition.
	i = dim - 1;
	r = dr * ((double)i - ghost + 0.5);
	// Fetch values.
	l_alpha = u[          i];
	l_psi   = u[    dim + i];
	phi     = u[2 * dim + i];
	Dr_l_alpha = Dr_u[          i];
	Dr_l_psi   = Dr_u[    dim + i];
	Dr_phi     = Dr_u[2 * dim + i];
	f[0 * dim + i] = -(r * Dr_l_alpha + l_alpha);
	f[1 * dim + i] = -(r * Dr_l_psi + l_psi);
	f[2 * dim + i] = -(r * Dr_phi + (1.0 + r * sqrt(m2 - w2)) * phi);

	// Finally set omega constraint RHS to zero.
	f[w_idx] = 0.0;

	// All done.
	return;
}
