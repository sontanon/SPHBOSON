#include "tools.h"
#include "param.h"

void diff1r(double *dvar, const double *var, const MKL_INT symr)
{
	// Auxiliar integer.
	MKL_INT i;

	// Inverse spatial step.
	double idr = 1.0 / dr;

	// Constants.
	const double half = 0.5;
	const double twelfth = 1.0 / 12.0;

	// Second-order derivatives.
	if (order == 2)
	{
		// Parity.
		dvar[ghost] = idr * half * (var[ghost + 1] - (double)symr * var[ghost]);
		#pragma omp parallel shared(dvar)
		{
			#pragma omp for schedule(guided)
			for (i = ghost + 1; i < dim - 1; i++)
			{
				dvar[i] = idr * half * (var[i + 1] - var[i - 1]);
			}
		}
		// Boundary
		dvar[dim - 1] = idr * half * (3.0 * var[dim - 1] - 4.0 * var[dim - 2] + var[dim - 3]);
	}
	// Fourth-order derivatives.
	else if (order == 4)
	{
		// Parity.
		dvar[ghost] = idr * twelfth * (-var[ghost + 2] + var[ghost + 1] * (8.0 + (double)symr) - (double)symr * 8.0 * var[ghost]);
		dvar[ghost + 1] = idr * twelfth * (-var[ghost + 3] + 8.0 * var[ghost + 2] + var[ghost] * (-8.0 + (double)symr));
		#pragma omp parallel shared(dvar)
		{
			#pragma omp for schedule(guided)
			for (i = ghost + 2; i < dim - 2; i++)
			{
				dvar[i] = idr * twelfth * (-var[i + 2] + 8.0 * var[i + 1] - 8.0 * var[i - 1] + var[i - 2]);
			}
		}
		// Boundary.
		dvar[dim - 2] = idr * twelfth * (3.0 * var[dim - 1] + 10.0 * var[dim - 2] - 18.0 * var[dim - 3] + 6.0 * var[dim - 4] - var[dim - 5]);
		dvar[dim - 1] = idr * twelfth * (25.0 * var[dim - 1] - 48.0 * var[dim - 2] + 36.0 * var[dim - 3] - 16.0 * var[dim - 4] + 3.0 * var[dim - 5]);
	}


	// Symmetries on axis.
	for (i = 0; i < ghost; i++)
	{
		dvar[ghost - 1 - i] = -(double)(symr) * dvar[ghost + i];
	}

	// All done.
	return;
}

void diff2r(double *dvar, const double *var, const MKL_INT symr)
{
	// Auxiliar integer.
	MKL_INT i;

	// Inverse spatial step.
	double idr2 = 1.0 / (dr * dr);

	// Constants.
	const double twelfth = 1.0 / 12.0;

	// Second-order derivatives.
	if (order == 2)
	{
		// Parity.
		dvar[ghost] = idr2 * (var[ghost + 1] + var[ghost] * (-2.0 + (double)symr));
		#pragma omp parallel shared(dvar)
		{
			#pragma omp for schedule(guided)
			for (i = ghost + 1; i < dim - 1; i++)
			{
				dvar[i] = idr2 * (var[i + 1] - 2.0 * var[i] + var[i - 1]);
			}
		}
		// Boundary.
		dvar[dim - 1] = idr2 * (2.0 * var[dim - 1] - 5.0 * var[dim - 2] + 4.0 * var[dim - 3] - var[dim - 4]);
	}
	// Fourth-order derivatives.
	else if (order == 4)
	{
		// Parity.
		dvar[ghost] = idr2 * twelfth * (-var[ghost + 2] + var[ghost + 1] * (16.0 - (double)symr) + var[ghost] * (-30.0 + 16.0 * (double)symr));
		dvar[ghost + 1] = idr2 * twelfth * (-var[ghost + 3] + 16.0 * var[ghost + 2] - 30.0 * var[ghost + 1] + var[ghost] * (16.0 - (double)symr)); 
		#pragma omp parallel shared(dvar)
		{
			#pragma omp for schedule(guided)
			for (i = ghost + 2; i < dim - 2; i++)
			{
				dvar[i] = idr2 * twelfth * (-var[i + 2] + 16.0 * var[i + 1] - 30.0 * var[i] + 16.0 * var[i - 1] - var[i - 2]);
			}
		}
		// Boundary.
		dvar[dim - 2] = idr2 * twelfth * (10.0 * var[dim - 1] - 15.0 * var[dim - 2] - 4.0 * var[dim - 3] + 14.0 * var[dim - 4] - 6.0 * var[dim - 5] + var[dim - 6]);
		dvar[dim - 1] = idr2 * twelfth * (45.0 * var[dim - 1] - 154.0 * var[dim - 2] + 214.0 * var[dim - 3] - 156.0 * var[dim - 4] + 61.0 * var[dim - 5] - 10.0 * var[dim - 6]);

	}

	// Symmetries on axis.
	for (i = 0; i < ghost; i++)
	{
		dvar[ghost - 1 - i] = (double)(symr) * dvar[ghost + i];
	}

	// All done.
	return;
}
