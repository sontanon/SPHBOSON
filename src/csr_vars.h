void jacobian_4th_order_variable_omega_c(
	double *aa,		// CSR array for values.
	MKL_INT *ia,		// CSR array for row beginnings.
	MKL_INT *ja,		// CSR array for columns.
	const MKL_INT dim,	// Grid total dimension.
	const MKL_INT ghost,	// Number of ghost zones.
	const MKL_INT i,	// Integer coordinate for r: 0 <= i < dim.
	const double dr,	// Spatial step for r. 
	const double m,		// Scalar field mass.
	const double xi,	// Scalar field frequency variable.
	// Now come the grid variables. For c stencil, each grid function has 5 variables.
	double u10, double u11, double u12, double u13, double u14,
	double u20, double u21, double u22, double u23, double u24,
	double u30, double u31, double u32, double u33, double u34,
	const MKL_INT offset1,
	const MKL_INT offset2,
	const MKL_INT offset3);

void jacobian_4th_order_variable_omega_s(
	double *aa,		// CSR array for values.
	MKL_INT *ia,		// CSR array for row beginnings.
	MKL_INT *ja,		// CSR array for columns.
	const MKL_INT dim,	// Grid total dimension.
	const MKL_INT ghost,	// Number of ghost zones.
	const MKL_INT i,	// Integer coordinate for r: 0 <= i < dim.
	const double dr,	// Spatial step for r. 
	const double m,		// Scalar field mass.
	const double xi,	// Scalar field frequency variable.
	// Now come the grid variables. For c stencil, each grid function has 6 variables.
	double u10, double u11, double u12, double u13, double u14, double u15,
	double u20, double u21, double u22, double u23, double u24, double u25,
	double u30, double u31, double u32, double u33, double u34, double u35,
	const MKL_INT offset1,
	const MKL_INT offset2,
	const MKL_INT offset3);