#include "tools.h"
#include "param.h"

double norm2(double *x)
{
	return cblas_dnrm2(dim, x, 1) / sqrt(dim);
}

double dot(double *x, double *y)
{
	return cblas_ddot(dim, x, 1, y, 1) / dim;
}
