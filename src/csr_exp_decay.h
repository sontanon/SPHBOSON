// Decay along r direction.
void exp_decay_4th_order(double *aa,
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
	const double m);