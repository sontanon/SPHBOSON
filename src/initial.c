#include "tools.h"
#include "param.h"
#include "omega_calc.h"

#define EPS 1.0

void initial_guess(double *u)
{
	// The main idea is to set:
	// 1. Lapse alpha to one => Logarithm to zero.
	// 2. Conformal factor to one => Logarithm to zero.
	// 3. Scalar field to gaussian configuration (see below).
	// 4. Omega to w0 => chi variable to inverse_omega_calc (see below).

	// Integer counter.
	MKL_INT i = 0, k = 0;

	// Auxiliary variables.
	double r;

	// Set omega variable.
	if (w_i)
	{
		read_single_file(&w0, w_i, 1, __FILE__, __LINE__);
		printf("***          Read omega initial data.       \n");
	}

	// Scale omega.
	w0 *= scale_u3;
	// Write chi variable.
	u[w_idx] = inverse_omega_calc(w0, m);

	//double m2 = m * m;
	//double w2 = w0 * w0;
	//double chi = sqrt(m2 - w2);

	if (readInitialData == 3)
	{
		double *u_0 = (double *)SAFE_MALLOC((GNUM * dim + 1) * sizeof(double));

		read_single_file(u_0 + 0 * dim, log_alpha_i	, dimInitial, __FILE__, __LINE__);
		read_single_file(u_0 + 1 * dim, log_psi_i	, dimInitial, __FILE__, __LINE__);
		read_single_file(u_0 + 2 * dim, phi_i		, dimInitial, __FILE__, __LINE__);

		u_0[GNUM * dimInitial] = u[w_idx];

		// Interpolate u0 into u.
		//initial_interpolator(u, u_0, dimInitial, ghost_i, order_i, dr_i, NrInterior, ghost, order, dr, w0, m);

		// Free initial data.
		SAFE_FREE(u_0);
	}
	else
	{
		if (!log_alpha_i)
		{
			#pragma omp parallel shared(u)
			{
				#pragma omp for schedule(guided)
				for (i = 0 * dim ; i < 1 * dim; i++)
					u[i] = 0.0;
			}
		}
		else
			read_single_file(u + 0 * dim, log_alpha_i, dim, __FILE__, __LINE__);

		if (!log_psi_i)
		{
			#pragma omp parallel shared(u)
			{
				#pragma omp for schedule(guided)
				for (i = 1 * dim ; i < 2 * dim; i++)
					u[i] = 0.0;
			}
		}
		else
			read_single_file(u + 1 * dim, log_psi_i, dim, __FILE__, __LINE__);

		if (!phi_i)
		{
			// Now do initial guess and calculate logarithm.
			// Initial guess is proportional to 
			//
			// phi = exp(-sqrt(m2 - w2) * r) * (EPS / (EPS + r)).
			//
			#pragma omp parallel shared(u) private(r)
			{
				#pragma omp for schedule(guided)
				for (i = 0; i < dim; i++)
				{
					r = dr * (i - ghost + 0.5);
					u[2* dim + i] = phi0 * exp(-0.5 * r * r / (sigma * sigma));
				}
			}
		}
		else
		{
			// Read file.
			read_single_file(u + 2 * dim, phi_i, dim, __FILE__, __LINE__);
		}
	}

	// Assert symmetries since they might not be automatic.
	for (k = 0; k < ghost; ++k)
	{
		u[0 * dim + k] = u[0 * dim + 2 * ghost - 1 - k];
		u[1 * dim + k] = u[1 * dim + 2 * ghost - 1 - k];
		u[2 * dim + k] = u[2 * dim + 2 * ghost - 1 - k];
	}

	// Before scaling, copy seed to u_seed.
	memcpy(u_seed + 0 * dim, u + 0 * dim, dim * sizeof(double));
	memcpy(u_seed + 1 * dim, u + 1 * dim, dim * sizeof(double));
	memcpy(u_seed + 2 * dim, u + 2 * dim, dim * sizeof(double));
	u_seed[w_idx] = u[w_idx];

	// Scale initial data.
	cblas_dscal(dim, scale_u0, u + 0 * dim, 1);
	cblas_dscal(dim, scale_u1, u + 1 * dim, 1);
	cblas_dscal(dim, scale_u2, u + 2 * dim, 1);
	// Omega has already been scaled.

	// All done.
	return;
}
