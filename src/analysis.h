void ex_phi_analysis(
	const MKL_INT print, 
	double *phi_max, 
	MKL_INT *f_res,
	double *u, 
	const MKL_INT ghost, 
	const MKL_INT order, 
	const MKL_INT dim);

void ex_analysis(const MKL_INT print,
	double *M_ADM,
	double *M_KOM,
	double *u,
	double *Dr_u,
	const MKL_INT ghost,
	const MKL_INT order,
	const MKL_INT dim, 
	const double dr);