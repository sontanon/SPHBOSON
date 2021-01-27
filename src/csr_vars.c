#include "tools.h"
#include "omega_calc.h"

#define GNUM 3

// Finite difference coefficients.
// Finite difference coefficients for 4th order.
#define D10 (+1.0 / 12.0)
#define D11 (-2.0 / 3.0)
#define D12 (+0.0)
#define D13 (+2.0 / 3.0)
#define D14 (-1.0 / 12.0)

#define D20 (-1.0 / 12.0)
#define D21 (+4.0 / 3.0)
#define D22 (-2.5)
#define D23 (+4.0 / 3.0)
#define D24 (-1.0 / 12.0)

#define S10 (+0.0)
#define S11 (-1.0 / 12.0)
#define S12 (+0.5)
#define S13 (-1.5)
#define S14 (+5.0 / 6.0)
#define S15 (+0.25)

#define S20 (+1.0 / 12.0)
#define S21 (-0.5)
#define S22 (+7.0 / 6.0)
#define S23 (-1.0 / 3.0)
#define S24 (-1.25)
#define S25 (+5.0 / 6.0)

// Jacobian for centered 4th order stencil and variable omega.
void jacobian_4th_order_variable_omega_c
(
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
	const MKL_INT offset3
)
{
	// Grid variables.
	double u1 = u12;
	double u2 = u22;
	double u3 = u32;

	// Physical names for readability.
	double alpha = exp(u1);
	double psi   = exp(u2);
	double phi   = u3;

	// Coordinates.
	double ri  = (double)i + 0.5 - ghost;
	double dr2 = dr * dr;

	// Scalar field mass and frequency.
	double w  = omega_calc(xi, m);
	double w2 = w * w;
	double m2 = m * m;

	// Omega variable index position.
	MKL_INT w_idx = GNUM * dim;

	// Short-hands.
	double phi2 = phi * phi;
	double alpha2 = alpha * alpha;
	double psi4 = psi * psi * psi * psi;

	// Finite differences.
	double dRu1 = D10 * u10 + D11 * u11 + D13 * u13 + D14 * u14;
	double dRu2 = D10 * u20 + D11 * u21 + D13 * u23 + D14 * u24;
	double dRu3 = D10 * u30 + D11 * u31 + D13 * u33 + D14 * u34;

	// Declare Jacobian submatrices.
	double jacobian_submatrix_1[3] = { 0.0 };  
	double jacobian_submatrix_2[3] = { 0.0 };  
	double jacobian_submatrix_3[3] = { 0.0 };  
	double jacobian_submatrix_w = 0.0;

	// CSR CODE FOR GRID NUMBER 1.
	
	// First write down Jacobian submatrices.
	// Submatrix 1.
	jacobian_submatrix_1[0] = 16.0 * M_PI * dr2 * psi4 * w2 * phi2 / alpha2;
	jacobian_submatrix_1[1] = 2.0 * (dRu1 + dRu2 + 1.0 / ri);
	jacobian_submatrix_1[2] = 1.0;

	// Submatrix 2.
	jacobian_submatrix_2[0] = 16.0 * M_PI * dr2 * psi4 * (m2 - 2.0 * w2 / alpha2) * phi2;
	jacobian_submatrix_2[1] = 2.0 * dRu1;
	jacobian_submatrix_2[2] = 0.0;

	// Submatrix 3.
	jacobian_submatrix_3[0] = 8.0 * M_PI * dr2 * psi4 * (m2 - 2.0 * w2 / alpha2) * phi;
	jacobian_submatrix_3[1] = 0.0;
	jacobian_submatrix_3[2] = 0.0;

	// Omega term.
	jacobian_submatrix_w = dw_du(xi, m) * (-16.0 * M_PI * dr2 * psi4 * w * phi2 / alpha2);

	// This row 0 * dim + i starts at offset1;
	ia[0 * dim + i] = BASE + offset1;

	// Values.
	aa[offset1 +  0] = D20*jacobian_submatrix_1[2]+D10*jacobian_submatrix_1[1];
	aa[offset1 +  1] = D21*jacobian_submatrix_1[2]+D11*jacobian_submatrix_1[1];
	aa[offset1 +  2] = D22*jacobian_submatrix_1[2]+1.0*jacobian_submatrix_1[0];
	aa[offset1 +  3] = D23*jacobian_submatrix_1[2]+D13*jacobian_submatrix_1[1];
	aa[offset1 +  4] = D24*jacobian_submatrix_1[2]+D14*jacobian_submatrix_1[1];

	aa[offset1 +  5] = D10*jacobian_submatrix_2[1];
	aa[offset1 +  6] = D11*jacobian_submatrix_2[1];
	aa[offset1 +  7] = 1.0*jacobian_submatrix_2[0];
	aa[offset1 +  8] = D13*jacobian_submatrix_2[1];
	aa[offset1 +  9] = D14*jacobian_submatrix_2[1];

	aa[offset1 + 10] = 1.0*jacobian_submatrix_3[0];

	aa[offset1 + 11] = jacobian_submatrix_w;

	// Columns.
	ja[offset1 +  0] = BASE + 0 * dim + i - 2;
	ja[offset1 +  1] = BASE + 0 * dim + i - 1;
	ja[offset1 +  2] = BASE + 0 * dim + i    ;
	ja[offset1 +  3] = BASE + 0 * dim + i + 1;
	ja[offset1 +  4] = BASE + 0 * dim + i + 2;

	ja[offset1 +  5] = BASE + 1 * dim + i - 2;
	ja[offset1 +  6] = BASE + 1 * dim + i - 1;
	ja[offset1 +  7] = BASE + 1 * dim + i    ;
	ja[offset1 +  8] = BASE + 1 * dim + i + 1;
	ja[offset1 +  9] = BASE + 1 * dim + i + 2;

	ja[offset1 + 10] = BASE + 2 * dim + i    ;

	ja[offset1 + 11] = BASE + w_idx;
	

	// CSR CODE FOR GRID NUMBER 2.
	
	// First write down Jacobian submatrices.
	// Submatrix 1.
	jacobian_submatrix_1[0] = -2.0 * M_PI * dr2 * psi4 * w2 * phi2 / alpha2;
	jacobian_submatrix_1[1] = 0.0;
	jacobian_submatrix_1[2] = 0.0;

	// Submatrix 2.
	jacobian_submatrix_2[0] = 4.0 * M_PI * dr2 * psi4 * (m2 + w2 / alpha2) * phi2;
	jacobian_submatrix_2[1] = 2.0 * (dRu2 + 1.0 / ri);
	jacobian_submatrix_2[2] = 1.0;

	// Submatrix 3.
	jacobian_submatrix_3[0] = 2.0 * M_PI * dr2 * psi4 * (m2 + w2 / alpha2) * phi;
	jacobian_submatrix_3[1] = 2.0 * M_PI * dRu3;
	jacobian_submatrix_3[2] = 0.0;

	// Omega term.
	jacobian_submatrix_w = dw_du(xi, m) * (2.0 * M_PI * dr2 * psi4 * w * phi2 / alpha2);

	// This row 1 * dim + i starts at offset2;
	ia[1 * dim + i] = BASE + offset2;

	// Values.
	aa[offset2 +  0] = 1.0*jacobian_submatrix_1[0];
	
	aa[offset2 +  1] = D20*jacobian_submatrix_2[2]+D10*jacobian_submatrix_2[1];
	aa[offset2 +  2] = D21*jacobian_submatrix_2[2]+D11*jacobian_submatrix_2[1];
	aa[offset2 +  3] = D22*jacobian_submatrix_2[2]+1.0*jacobian_submatrix_2[0];
	aa[offset2 +  4] = D23*jacobian_submatrix_2[2]+D13*jacobian_submatrix_2[1];
	aa[offset2 +  5] = D24*jacobian_submatrix_2[2]+D14*jacobian_submatrix_2[1];

	aa[offset2 +  6] = D10*jacobian_submatrix_3[1];
	aa[offset2 +  7] = D11*jacobian_submatrix_3[1];
	aa[offset2 +  8] = 1.0*jacobian_submatrix_3[0];
	aa[offset2 +  9] = D13*jacobian_submatrix_3[1];
	aa[offset2 + 10] = D14*jacobian_submatrix_3[1];

	aa[offset2 + 11] = jacobian_submatrix_w;

	// Columns.
	ja[offset2 +  0] = BASE + 0 * dim + i    ;

	ja[offset2 +  1] = BASE + 1 * dim + i - 2;
	ja[offset2 +  2] = BASE + 1 * dim + i - 1;
	ja[offset2 +  3] = BASE + 1 * dim + i    ;
	ja[offset2 +  4] = BASE + 1 * dim + i + 1;
	ja[offset2 +  5] = BASE + 1 * dim + i + 2;

	ja[offset2 +  6] = BASE + 2 * dim + i - 2;
	ja[offset2 +  7] = BASE + 2 * dim + i - 1;
	ja[offset2 +  8] = BASE + 2 * dim + i    ;
	ja[offset2 +  9] = BASE + 2 * dim + i + 1;
	ja[offset2 + 10] = BASE + 2 * dim + i + 2;

	ja[offset2 + 11] = BASE + w_idx;


	// CSR CODE FOR GRID NUMBER 3.
	
	// First write down Jacobian submatrices.
	// Submatrix 1.
	jacobian_submatrix_1[0] = -2.0 * dr2 * psi4 * w2 * phi / alpha2;
	jacobian_submatrix_1[1] = dRu3;
	jacobian_submatrix_1[2] = 0.0;

	// Submatrix 2.
	jacobian_submatrix_2[0] = 4.0 * dr2 * psi4 * (w2 / alpha2 - m2) * phi;
	jacobian_submatrix_2[1] = 2.0 * dRu3;
	jacobian_submatrix_2[2] = 0.0;

	// Submatrix 3.
	jacobian_submatrix_3[0] = dr2 * psi4 * (w2 / alpha - m2);
	jacobian_submatrix_3[1] = dRu1 + 2.0 * (dRu2 + 1.0 / ri);
	jacobian_submatrix_3[2] = 1.0;

	// Omega term.
	jacobian_submatrix_w = dw_du(xi, m) * (2.0 * dr2 * psi4 * w * phi / alpha2);

	// This row 2 * dim + i starts at offset3;
	ia[2 * dim + i] = BASE + offset3;

	// Values.
	aa[offset3 +  0] = D10*jacobian_submatrix_1[1];
	aa[offset3 +  1] = D11*jacobian_submatrix_1[1];
	aa[offset3 +  2] = 1.0*jacobian_submatrix_1[0];
	aa[offset3 +  3] = D13*jacobian_submatrix_1[1];
	aa[offset3 +  4] = D14*jacobian_submatrix_1[1];

	aa[offset3 +  5] = D10*jacobian_submatrix_2[1];
	aa[offset3 +  6] = D11*jacobian_submatrix_2[1];
	aa[offset3 +  7] = 1.0*jacobian_submatrix_2[0];
	aa[offset3 +  8] = D13*jacobian_submatrix_2[1];
	aa[offset3 +  9] = D14*jacobian_submatrix_2[1];

	aa[offset3 + 10] = D20*jacobian_submatrix_3[2]+D10*jacobian_submatrix_3[1];
	aa[offset3 + 11] = D21*jacobian_submatrix_3[2]+D11*jacobian_submatrix_3[1];
	aa[offset3 + 12] = D22*jacobian_submatrix_3[2]+1.0*jacobian_submatrix_3[0];
	aa[offset3 + 13] = D23*jacobian_submatrix_3[2]+D13*jacobian_submatrix_3[1];
	aa[offset3 + 14] = D24*jacobian_submatrix_3[2]+D14*jacobian_submatrix_3[1];

	aa[offset3 + 15] = jacobian_submatrix_w;

	// Columns.
	ja[offset3 +  0] = BASE + 0 * dim + i - 2;
	ja[offset3 +  1] = BASE + 0 * dim + i - 1;
	ja[offset3 +  2] = BASE + 0 * dim + i    ;
	ja[offset3 +  3] = BASE + 0 * dim + i + 1;
	ja[offset3 +  4] = BASE + 0 * dim + i + 2;

	ja[offset3 +  5] = BASE + 1 * dim + i - 2;
	ja[offset3 +  6] = BASE + 1 * dim + i - 1;
	ja[offset3 +  7] = BASE + 1 * dim + i    ;
	ja[offset3 +  8] = BASE + 1 * dim + i + 1;
	ja[offset3 +  9] = BASE + 1 * dim + i + 2;

	ja[offset3 + 10] = BASE + 2 * dim + i - 2;
	ja[offset3 + 11] = BASE + 2 * dim + i - 1;
	ja[offset3 + 12] = BASE + 2 * dim + i    ;
	ja[offset3 + 13] = BASE + 2 * dim + i + 1;
	ja[offset3 + 14] = BASE + 2 * dim + i + 2;

	ja[offset3 + 15] = BASE + w_idx;


	// All done.
	return;
}

// Jacobian for semi-onesided 4th order stencil and variable omega.
void jacobian_4th_order_variable_omega_s
(
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
	const MKL_INT offset3
)
{
	// Grid variables.
	double u1 = u14;
	double u2 = u24;
	double u3 = u34;

	// Physical names for readability.
	double alpha = exp(u1);
	double psi   = exp(u2);
	double phi   = u3;

	// Coordinates.
	double ri  = (double)i + 0.5 - ghost;
	//double r   = ri * dr;
	double dr2 = dr * dr;

	// Scalar field mass and frequency.
	double w  = omega_calc(xi, m);
	double w2 = w * w;
	double m2 = m * m;

	// Omega variable index position.
	MKL_INT w_idx = GNUM * dim;

	// Short-hands.
	double phi2 = phi * phi;
	double alpha2 = alpha * alpha;
	double psi4 = psi * psi * psi * psi;

	// Finite differences.
	double dRu1 = S11 * u11 + S12 * u12 + S13 * u13 + S14 * u14 + S15 * u15;
	double dRu2 = S11 * u21 + S12 * u22 + S13 * u23 + S14 * u24 + S15 * u25;
	double dRu3 = S11 * u31 + S12 * u32 + S13 * u33 + S14 * u34 + S15 * u35;

	// Declare Jacobian submatrices.
	double jacobian_submatrix_1[3] = { 0.0 };  
	double jacobian_submatrix_2[3] = { 0.0 };  
	double jacobian_submatrix_3[3] = { 0.0 };  
	double jacobian_submatrix_w = 0.0;

	// CSR CODE FOR GRID NUMBER 1.
	
	// First write down Jacobian submatrices.
	// Submatrix 1.
	jacobian_submatrix_1[0] = 16.0 * M_PI * dr2 * psi4 * w2 * phi2 / alpha2;
	jacobian_submatrix_1[1] = 2.0 * (dRu1 + dRu2 + 1.0 / ri);
	jacobian_submatrix_1[2] = 1.0;

	// Submatrix 2.
	jacobian_submatrix_2[0] = 16.0 * M_PI * dr2 * psi4 * (m2 - 2.0 * w2 / alpha2) * phi2;
	jacobian_submatrix_2[1] = 2.0 * dRu1;
	jacobian_submatrix_2[2] = 0.0;

	// Submatrix 3.
	jacobian_submatrix_3[0] = 8.0 * M_PI * dr2 * psi4 * (m2 - 2.0 * w2 / alpha2) * phi;
	jacobian_submatrix_3[1] = 0.0;
	jacobian_submatrix_3[2] = 0.0;

	// Omega term.
	jacobian_submatrix_w = dw_du(xi, m) * (-16.0 * M_PI * dr2 * psi4 * w * phi2 / alpha2);

	// This row 0 * dim + i starts at offset1;
	ia[0 * dim + i] = BASE + offset1;

	// Values.
	aa[offset1 +  0] = S20*jacobian_submatrix_1[2];
	aa[offset1 +  1] = S21*jacobian_submatrix_1[2]+S11*jacobian_submatrix_1[1];
	aa[offset1 +  2] = S22*jacobian_submatrix_1[2]+S12*jacobian_submatrix_1[1];
	aa[offset1 +  3] = S23*jacobian_submatrix_1[2]+S13*jacobian_submatrix_1[1];
	aa[offset1 +  4] = S24*jacobian_submatrix_1[2]+S14*jacobian_submatrix_1[1]+1.0*jacobian_submatrix_1[0];
	aa[offset1 +  5] = S25*jacobian_submatrix_1[2]+S15*jacobian_submatrix_1[1];

	aa[offset1 +  6] = S11*jacobian_submatrix_2[1];
	aa[offset1 +  7] = S12*jacobian_submatrix_2[1];
	aa[offset1 +  8] = S13*jacobian_submatrix_2[1];
	aa[offset1 +  9] = S14*jacobian_submatrix_2[1]+1.0*jacobian_submatrix_2[0];
	aa[offset1 + 10] = S15*jacobian_submatrix_2[1];

	aa[offset1 + 11] = 1.0*jacobian_submatrix_3[0];

	aa[offset1 + 12] = jacobian_submatrix_w;

	// Columns.
	ja[offset1 +  0] = BASE + 0 * dim + i - 4;
	ja[offset1 +  1] = BASE + 0 * dim + i - 3;
	ja[offset1 +  2] = BASE + 0 * dim + i - 2;
	ja[offset1 +  3] = BASE + 0 * dim + i - 1;
	ja[offset1 +  4] = BASE + 0 * dim + i    ;
	ja[offset1 +  5] = BASE + 0 * dim + i + 1;

	ja[offset1 +  6] = BASE + 1 * dim + i - 3;
	ja[offset1 +  7] = BASE + 1 * dim + i - 2;
	ja[offset1 +  8] = BASE + 1 * dim + i - 1;
	ja[offset1 +  9] = BASE + 1 * dim + i    ;
	ja[offset1 + 10] = BASE + 1 * dim + i + 1;

	ja[offset1 + 11] = BASE + 2 * dim + i    ;

	ja[offset1 + 12] = BASE + w_idx;
	

	// CSR CODE FOR GRID NUMBER 2.
	
	// First write down Jacobian submatrices.
	// Submatrix 1.
	jacobian_submatrix_1[0] = -2.0 * M_PI * dr2 * psi4 * w2 * phi2 / alpha2;
	jacobian_submatrix_1[1] = 0.0;
	jacobian_submatrix_1[2] = 0.0;

	// Submatrix 2.
	jacobian_submatrix_2[0] = 4.0 * M_PI * dr2 * psi4 * (m2 + w2 / alpha2) * phi2;
	jacobian_submatrix_2[1] = 2.0 * (dRu2 + 1.0 / ri);
	jacobian_submatrix_2[2] = 1.0;

	// Submatrix 3.
	jacobian_submatrix_3[0] = 2.0 * M_PI * dr2 * psi4 * (m2 + w2 / alpha2) * phi;
	jacobian_submatrix_3[1] = 2.0 * M_PI * dRu3;
	jacobian_submatrix_3[2] = 0.0;

	// Omega term.
	jacobian_submatrix_w = dw_du(xi, m) * (2.0 * M_PI * dr2 * psi4 * w * phi2 / alpha2);

	// This row 1 * dim + i starts at offset2;
	ia[1 * dim + i] = BASE + offset2;

	// Values.
	aa[offset2 +  0] = 1.0*jacobian_submatrix_1[0];

	aa[offset2 +  1] = S20*jacobian_submatrix_2[2];
	aa[offset2 +  2] = S21*jacobian_submatrix_2[2]+S11*jacobian_submatrix_2[1];
	aa[offset2 +  3] = S22*jacobian_submatrix_2[2]+S12*jacobian_submatrix_2[1];
	aa[offset2 +  4] = S23*jacobian_submatrix_2[2]+S13*jacobian_submatrix_2[1];
	aa[offset2 +  5] = S24*jacobian_submatrix_2[2]+S14*jacobian_submatrix_2[1]+1.0*jacobian_submatrix_2[0];
	aa[offset2 +  6] = S25*jacobian_submatrix_2[2]+S15*jacobian_submatrix_2[1];
	
	aa[offset2 +  7] = S11*jacobian_submatrix_3[1];
	aa[offset2 +  8] = S12*jacobian_submatrix_3[1];
	aa[offset2 +  9] = S13*jacobian_submatrix_3[1];
	aa[offset2 + 10] = S14*jacobian_submatrix_3[1]+1.0*jacobian_submatrix_3[0];
	aa[offset2 + 11] = S15*jacobian_submatrix_3[1];

	aa[offset2 + 12] = jacobian_submatrix_w;

	// Columns.
	ja[offset2 +  0] = BASE + 0 * dim + i    ;

	ja[offset2 +  1] = BASE + 1 * dim + i - 4;
	ja[offset2 +  2] = BASE + 1 * dim + i - 3;
	ja[offset2 +  3] = BASE + 1 * dim + i - 2;
	ja[offset2 +  4] = BASE + 1 * dim + i - 1;
	ja[offset2 +  5] = BASE + 1 * dim + i    ;
	ja[offset2 +  6] = BASE + 1 * dim + i + 1;

	ja[offset2 +  7] = BASE + 2 * dim + i - 3;
	ja[offset2 +  8] = BASE + 2 * dim + i - 2;
	ja[offset2 +  9] = BASE + 2 * dim + i - 1;
	ja[offset2 + 10] = BASE + 2 * dim + i    ;
	ja[offset2 + 11] = BASE + 2 * dim + i + 1;

	ja[offset2 + 12] = BASE + w_idx;


	// CSR CODE FOR GRID NUMBER 3.
	
	// First write down Jacobian submatrices.
	// Submatrix 1.
	jacobian_submatrix_1[0] = -2.0 * dr2 * psi4 * w2 * phi / alpha2;
	jacobian_submatrix_1[1] = dRu3;
	jacobian_submatrix_1[2] = 0.0;

	// Submatrix 2.
	jacobian_submatrix_2[0] = 4.0 * dr2 * psi4 * (w2 / alpha2 - m2) * phi;
	jacobian_submatrix_2[1] = 2.0 * dRu3;
	jacobian_submatrix_2[2] = 0.0;

	// Submatrix 3.
	jacobian_submatrix_3[0] = dr2 * psi4 * (w2 / alpha - m2);
	jacobian_submatrix_3[1] = dRu1 + 2.0 * (dRu2 + 1.0 / ri);
	jacobian_submatrix_3[2] = 1.0;

	// Omega term.
	jacobian_submatrix_w = dw_du(xi, m) * (2.0 * dr2 * psi4 * w * phi / alpha2);

	// This row 2 * dim + i starts at offset3;
	ia[2 * dim + i] = BASE + offset3;

	// Values.
	aa[offset3 +  0] = S11*jacobian_submatrix_1[1];
	aa[offset3 +  1] = S12*jacobian_submatrix_1[1];
	aa[offset3 +  2] = S13*jacobian_submatrix_1[1];
	aa[offset3 +  3] = S14*jacobian_submatrix_1[1]+1.0*jacobian_submatrix_1[0];
	aa[offset3 +  4] = S15*jacobian_submatrix_1[1];

	aa[offset3 +  5] = S11*jacobian_submatrix_2[1];
	aa[offset3 +  6] = S12*jacobian_submatrix_2[1];
	aa[offset3 +  7] = S13*jacobian_submatrix_2[1];
	aa[offset3 +  8] = S14*jacobian_submatrix_2[1]+1.0*jacobian_submatrix_2[0];
	aa[offset3 +  9] = S15*jacobian_submatrix_2[1];

	aa[offset3 + 10] = S20*jacobian_submatrix_3[2];
	aa[offset3 + 11] = S21*jacobian_submatrix_3[2]+S11*jacobian_submatrix_3[1];
	aa[offset3 + 12] = S22*jacobian_submatrix_3[2]+S12*jacobian_submatrix_3[1];
	aa[offset3 + 13] = S23*jacobian_submatrix_3[2]+S13*jacobian_submatrix_3[1];
	aa[offset3 + 14] = S24*jacobian_submatrix_3[2]+S14*jacobian_submatrix_3[1]+1.0*jacobian_submatrix_3[0];
	aa[offset3 + 15] = S25*jacobian_submatrix_3[2]+S15*jacobian_submatrix_3[1];

	aa[offset3 + 16] = jacobian_submatrix_w;

	// Columns.
	ja[offset3 +  0] = BASE + 0 * dim + i - 3;
	ja[offset3 +  1] = BASE + 0 * dim + i - 2;
	ja[offset3 +  2] = BASE + 0 * dim + i - 1;
	ja[offset3 +  3] = BASE + 0 * dim + i    ;
	ja[offset3 +  4] = BASE + 0 * dim + i + 1;

	ja[offset3 +  5] = BASE + 1 * dim + i - 3;
	ja[offset3 +  6] = BASE + 1 * dim + i - 2;
	ja[offset3 +  7] = BASE + 1 * dim + i - 1;
	ja[offset3 +  8] = BASE + 1 * dim + i    ;
	ja[offset3 +  9] = BASE + 1 * dim + i + 1;

	ja[offset3 + 10] = BASE + 2 * dim + i - 4;
	ja[offset3 + 11] = BASE + 2 * dim + i - 3;
	ja[offset3 + 12] = BASE + 2 * dim + i - 2;
	ja[offset3 + 13] = BASE + 2 * dim + i - 1;
	ja[offset3 + 14] = BASE + 2 * dim + i    ;
	ja[offset3 + 15] = BASE + 2 * dim + i + 1;

	ja[offset3 + 16] = BASE + w_idx;


	// All done.
	return;
}