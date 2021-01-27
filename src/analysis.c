#include "tools.h"
#include "derivatives.h"
#include "simpson.h"

#define EVEN 1

void ex_phi_analysis(
	const MKL_INT print, 
	double *phi_max, 
	MKL_INT *f_res,
	double *u, 
	const MKL_INT ghost, 
	const MKL_INT order, 
	const MKL_INT dim)
{
	MKL_INT k = 0;

	// Extract variable.
	double *phi = u + 2 * dim;

	double phi_min = phi[cblas_idamin(dim, phi, 1)];

	// Parabolic interpolation.
	*phi_max = (41.0 / 48.0) * phi[ghost] + (5.0 / 32.0) * phi[ghost + 1] - (1.0 / 96.0) * phi[ghost + 2];

	// Count number of points at half-width length.
	double hwl = 0.5 * *phi_max;
	*f_res = 0;

	for (k = 1; k < dim; ++k)
	{
		if (phi[k] > hwl)
			++(*f_res);
		else
			break;
	}

	if (print)
	{
		write_single_file(phi_max, "phi_max.asc", 1);
		write_single_file_integer(f_res, "hwl_resolution.asc", 1);
		printf("***\n");
		printf("*** Scalar Field Analysis\n");
		printf("***\n");
		printf("***  -------------------------- ----------------------- ---------------- \n");
		printf("*** | max(phi) or phi(0)       | min(phi)              | HWL Resolution |\n");
		printf("***  -------------------------- ----------------------- ---------------- \n");	
		printf("*** |       %-6.5e        |      %-6.5e      |      % 4lld      |\n", *phi_max, phi_min, *f_res);
		printf("***  -------------------------- ----------------------- ---------------- \n");
		printf("***\n");
	}

	// All done.
	return;
}

void ex_analysis(const MKL_INT print,
	double *M_ADM,
	double *M_KOM,
	double *u,
	double *Dr_u,
	const MKL_INT ghost,
	const MKL_INT order,
	const MKL_INT dim, 
	const double dr)
{
	// Masses.
	double *array_M_ADM = (double *)SAFE_MALLOC(sizeof(double) * dim);
	double *array_M_KOM = (double *)SAFE_MALLOC(sizeof(double) * dim);

	MKL_INT i = 0;
	double psi, psi2, psi4;
	double r, r2;

	diff1r(Dr_u + 0 * dim, u + 0 * dim, 1);
	diff1r(Dr_u + 1 * dim, u + 1 * dim, 1);

	for (i = 0; i < dim; ++i)
	{
		r = dr * (i - ghost + 0.5);
		r2 = r * r;
		psi = exp(u[1 * dim + i]);
		psi2 = psi * psi;
		psi4 = psi2 * psi2;
		// M = -2.0 * psi4 * r2 * dR(log(psi)).
		array_M_ADM[i] = -2.0 * psi4 * r2 * Dr_u[1 * dim + i];
		// M = psi2 * r2 * dR(log(alpha)) * alpha.
		array_M_KOM[i] = r2 * psi2 * exp(u[0 * dim + i]) * Dr_u[0 * dim + i];
	}

	*M_ADM = array_M_ADM[dim - 1];
	*M_KOM = array_M_KOM[dim - 1];

	if (print)
	{
		// Print information to screen.
		printf("*** \n");
		printf("*** GLOBAL QUANTITIES ANALYSIS\n");
		printf("*** \n");
		printf("*** Final radius is r_inf = %6.5e.\n", dr * (dim - ghost - 0.5));
		printf("*** \n");
		printf("***  ----------------------- ----------------- \n");
		printf("*** | Komar Mass            | ADM Mass        |\n");
		printf("***  ----------------------- ----------------- \n");	
		//printf("*** |      1234567890123       |     1234567890123     |  1234567890123  |  1234567890123  |      1234567890123       |     1234567890123      |\n");
		printf("*** |      %-6.5e      |   %-6.5e   |\n", *M_KOM, *M_ADM);
		printf("***  ----------------------- ----------------- \n");
		printf("*** \n");

		write_single_file(array_M_ADM, "M_ADM.asc", dim);
		write_single_file(array_M_KOM, "M_Komar.asc", dim);
	}

	SAFE_FREE(array_M_ADM);
	SAFE_FREE(array_M_KOM);

	// All done.
	return;
}