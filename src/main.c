#include "tools.h"
#define MAIN_FILE
#include "param.h"

#include "parser.h"
#include "io.h"
#include "initial.h"
#include "rhs.h"
#include "pardiso_start.h"
#include "pardiso_stop.h"
#include "omega_calc.h"
#include "csr.h"
#include "nleq_err.h"
#include "nleq_res.h"
#include "vector_algebra.h"
#include "pardiso_solve.h"
#include "analysis.h"
#include "low_rank.h"
#include "newton.h"

int main(int argc, char *argv[])
{
	// I/O error code.
	int io_code = 0;

	// Integer counter.
	MKL_INT i = 0, k = 0;

	MKL_INT counter_i = 0;

	// Error code.
	MKL_INT errCode = 1;

	// Initial message.
	printf("***************************************************\n");
	printf("***************************************************\n");
	printf("***                                                \n");
	printf("***                    SPHBOSON                    \n");
	printf("***                                                \n");
	printf("***          Global Newton Method Version          \n");
	printf("***                                                \n");
	printf("***        Author: Santiago Ontanon Sanchez        \n");
	printf("***                                                \n");
	printf("***              ICN UNAM, Mexico City             \n");
	printf("***                                                \n");
	printf("***                                                \n");
	printf("***             First Revision: 8/2019             \n");
	printf("***                                                \n");
	printf("***             Last  Revision: 1/2021             \n");
	printf("***                                                \n");
	printf("***************************************************\n");

	// File name is in argv[1]. Check that we have at
	// least one argument.
	if (argc < 2)
	{
		printf("***                                                \n");
		printf("***           Usage: ./SPHBOSON file.par           \n");
		printf("***                                                \n");
		printf("***            Missing parameter file.             \n");
		printf("***                                                \n");
		printf("***************************************************\n");
		printf("***************************************************\n");
		return EXIT_FAILURE;
	}

	// Parse arguments.
	parser(argv[1]);

	// Future scale factors.
	double next_scale[GNUM + 1] = { scale_u0, scale_u1, scale_u2, scale_u3 };

	// Peaks.
	double peak_next[GNUM + 1] = { 0.0 };
	double peak_prev[GNUM + 1] = { 0.0 };

	// Print main arguements.
	printf("***************************************************\n");
	printf("***                                                \n");
	printf("***           Generating Spherical Boson           \n");
	printf("***            Star Initial Data For NR.           \n");
	printf("***                                                \n");
	printf("***           GRID:                                \n");
	printf("***            dr          = %-12.10lf          \n", dr);
	printf("***            dim         = %-7lld               \n", dim);
	printf("***            NrInterior  = %-7lld               \n", NrInterior);
	printf("***            order       = %lld                     \n", order);
	printf("***            ghost       = %lld                     \n", ghost);
	printf("***                                                \n");
	printf("***           SCALAR FIELD:                        \n");
	printf("***            m           = %-12.10lf          \n", m);
	printf("***            phi0        = %-12.10lf          \n", phi0);
	printf("***            w0          = %-12.10lf          \n", w0);
	if (fixedPhi)
	{
		printf("***            r(fixedPhi) = %-12.10lf          \n", dr * (fixedPhi - 0.5));
	}
	else
	{
		printf("***            Initial Omega is Fixed.             \n");
	}
	printf("***                                                \n");
	printf("***           INITIAL DATA:                        \n");
	printf("***            readInitialData = %lld     \n", readInitialData);
	if (readInitialData)
	{
		printf("***            log_alpha_i     = %-18s     \n", log_alpha_i);
		printf("***            log_psi_i       = %-18s     \n", log_psi_i);
		printf("***            phi_i           = %-18s     \n", phi_i);
		printf("***            w_i             = %-18s     \n", w_i);
	}
	else
	{
		printf("***            phi0            = %-12.10E  \n", phi0);
	}
	if (!w_i)
	{
		printf("***            w0          = %-12.10E          \n", w0);
	}
	printf("***                                                \n");
	printf("***           SOLVER:                              \n");
	printf("***            solverType    = %-18s  \n", (solverType == 1) ? "Error" : "Residual");
	printf("***            epsilon       = %-12.10E        \n", epsilon);
	printf("***            maxNewtonIter = %-4lld                \n", maxNewtonIter);
	printf("***            lambda0       = %-12.10E        \n", lambda0);
	printf("***            lambdaMin     = %-12.10E        \n", lambdaMin);
	printf("***            useLowRank    = %lld       \n", useLowRank);
	printf("***                                                \n");
	printf("******************************************************\n");

	#pragma omp parallel
	{
		#pragma omp master
		{
			// Determine OMP threads.
			printf("***************************************************\n");
			printf("***                                                \n");
			printf("***            Maximum OMP threads = %d             \n", omp_get_max_threads());
			printf("***            Currently running on %d              \n", omp_get_num_threads());
			printf("***                                                \n");
			printf("***************************************************\n");
		}
	}


	printf("***************************************************\n");
	printf("***                                                \n");
	printf("***               Allocating memory...             \n");
	printf("***                                                \n");


	// Allocate pointer to double pointers.
	double **u      = (double **)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double *));
	double **f      = (double **)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double *));
	double **du     = (double **)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double *));
	double **du_bar = (double **)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double *));

	// Allocate memory.
	for (i = 0; i < maxNewtonIter; i++)
	{
		u[i]      = (double *)SAFE_MALLOC((GNUM * dim + 1) * sizeof(double));
		f[i]      = (double *)SAFE_MALLOC((GNUM * dim + 1) * sizeof(double));
		du[i]     = (double *)SAFE_MALLOC((GNUM * dim + 1) * sizeof(double));
		du_bar[i] = (double *)SAFE_MALLOC((GNUM * dim + 1) * sizeof(double));
	}

	// Also include radial grid.
	double *r = (double *)SAFE_MALLOC(dim * sizeof(double));

	// Initial data seed.
	u_seed = (double *)SAFE_MALLOC(dim * sizeof(double));

	// Auxiliary derivative integers.
	Dr_u  = (double *)SAFE_MALLOC((3 * dim + 1) * sizeof(double));
	Drr_u = (double *)SAFE_MALLOC((3 * dim + 1) * sizeof(double));

	// Newton output parameters.
	double *norm_f		= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));
	double *norm_du		= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));
	double *norm_du_bar	= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));
	double *lambda	= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));
	double *Theta	= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));
	double *mu	= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));
	double *lambda_prime	= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));
	double *mu_prime	= (double *)SAFE_MALLOC((maxNewtonIter + 1) * sizeof(double));

	// Initial guess norms.
	double f_norms[GNUM];

	// Final omega.
	double w = m;

	// Fill r grid.
	#pragma omp parallel shared(r)
	{
		#pragma omp for schedule(guided)
		for (i = 0; i < dim; i++)
		{
			r[i] = ((double)(i - ghost) + 0.5) * dr;
		}
	}

	printf("***               Finished allocation!             \n");
	printf("***                                                \n");
	printf("***************************************************\n");

	printf("***************************************************\n");
	printf("***                                                \n");
	printf("***           Allocating PARDISO memory...         \n");
	printf("***                                                \n");

	// Initialize PARDISO memory and paramters.
	// Square matrix dimension is (GNUM * dim + 1).
	pardiso_start(GNUM * dim + 1);

	// Allocate CSR matrix.
	csr_matrix J;
	MKL_INT nnz = nnz_jacobian();
	csr_allocate(&J, GNUM * dim + 1, GNUM * dim + 1, nnz);

	printf("***                                                \n");
	printf("***            Allocated CSR matrix with:          \n");
	printf("***             Rows      = %-6lld                 \n", J.nrows);
	printf("***             Columns   = %-6lld                 \n", J.ncols);
	printf("***             Non-zeros = %-12lld           \n", J.nnz);
	printf("***                                                \n");
	printf("***           Finished PARDISO allocation!         \n");
	printf("***                                                \n");
	printf("***************************************************\n");

	// LOW RANK UPDATE and linear solver subroutines.
	// Linear Solver Subroutine.
	void (*linear_solve_1)(double *, csr_matrix *, double *);
	if (useLowRank)
	{
		linear_solve_1 = pardiso_solve_low_rank;
		diff_gen();
	}
	else
		linear_solve_1 = pardiso_simple_solve;

	/* NO LOW-RANK UPDATE
	void (*linear_solve_1)(double *, csr_matrix *, double *);
	linear_solve_1 = pardiso_simple_solve;
	*/
	
	// Repeated solver.
	void (*linear_solve_2)(double *, csr_matrix *, double *);
	linear_solve_2 = pardiso_repeated_solve;

	printf("***************************************************\n");
	printf("***                                                \n");
	printf("***          Setting initial guess and RHS.        \n");

	// Set initial guess.
	initial_guess(u[0]);

	// Loop over sweep.
	while (1)
	{
		// DO I/O.
		io(initial_dirname, argv[1]);

		// Print main variables.
		write_single_file(u[0]          , "log_alpha_i.asc", dim);
		write_single_file(u[0] +     dim, "log_psi_i.asc", dim);
		write_single_file(u[0] + 2 * dim, "phi_i.asc", dim);
		write_single_file(&w0, "w_i.asc", 1);

		// Write coordinate grids and print where phi is fixed.
		write_single_file(r, "r.asc", dim);

		// And initial "seed."
		write_single_file(u_seed          , "log_alpha_seed.asc", dim);
		write_single_file(u_seed +     dim, "log_psi_seed.asc", dim);
		write_single_file(u_seed + 2 * dim, "phi_seed.asc", dim);
		w = omega_calc(u_seed[GNUM * dim], m);
		write_single_file(&w, "w_seed.asc", 1);

		// First calculate initial RHS.
		rhs(f[0], u[0]);

		// Print initial RHS.
		write_single_file(f[0]          , "f0_i.asc", dim);
		write_single_file(f[0] +     dim, "f1_i.asc", dim);
		write_single_file(f[0] + 2 * dim, "f2_i.asc", dim);

		// Calculate 2-norms.
		f_norms[0] = norm2(f[0]          );
		f_norms[1] = norm2(f[0] +     dim);
		f_norms[2] = norm2(f[0] + 2 * dim);
		
		printf("***                                                \n");
		printf("***        INITIAL GUESS:                          \n");
		printf("***           || f0 ||   = %-12.10E           \n", f_norms[0]);
		printf("***           || f1 ||   = %-12.10E           \n", f_norms[1]);
		printf("***           || f2 ||   = %-12.10E           \n", f_norms[2]);
		printf("***                                                \n");
		printf("******************************************************\n");

		// Set initial damping factor lambda[0].
		lambda[0] = lambda0;

		// Set k = 0.
		k = 0;

		/* MAIN ALGORITHM: NEWTON SOLVER */
		// Start Newton iterations.
		if (maxNewtonIter > 0)
		{
			switch (solverType)
			{
				// Error-based algorithm.
				case 1:
					// Call algorithm.
					k = nleq_err(&errCode, u, f, lambda,
							du, du_bar, norm_du, norm_du_bar,
							Theta, mu, lambda_prime, mu_prime,
							&J, epsilon, maxNewtonIter, 8, 8,
							lambdaMin, localSolver,
							rhs, csr_gen_jacobian, 
							norm2, dot,
							linear_solve_1, linear_solve_2);
					break;
				// Residual-based algorithm.
				case 2:
					// ||f[0]|| is also an input parameter.
					norm_f[0] = norm2(u[0]);
					// Calle algorithm.
					k = nleq_res(&errCode, u, f, lambda,
							du, norm_f, Theta, mu, lambda_prime, mu_prime,
							&J, epsilon, maxNewtonIter, 8, 8,
							lambdaMin, localSolver, 
							rhs, csr_gen_jacobian, 
							norm2, dot,
							linear_solve_1, linear_solve_2);
					break;
				case 3: 
					// ||f[0]|| is also an input parameter.
					norm_f[0] = norm2(u[0]);
					k = newton(&errCode, u, f, lambda,
							du, norm_du, Theta,
							&J, epsilon, maxNewtonIter,
							rhs, csr_gen_jacobian,
							norm2,
							linear_solve_1);
					break;
			}

			// Write errCode to file.
			write_single_file_integer(&errCode, "error_code.asc", 1);

			// Check for convergence.
			if (errCode != 0)
			{
				printf("******************************************************\n");
				printf("***                                                \n");
				printf("***    Warning! Did not converge: Error Code = %lld  \n", errCode);
				printf("***    Will output anyway. Do not trust results!   \n");
				printf("***                                                \n");
				printf("******************************************************\n");
				k = -k;
			}
		}
		else
		{
			printf("******************************************************\n");
			printf("***                                                \n");
			printf("***    Warning! User did not specify any Newton Iterations.  \n");
			printf("***                                                \n");
			printf("******************************************************\n");
			k = 0;
		}

		// Get omega.
		w = omega_calc(u[k][w_idx], m);

		// Print final solutions
		write_single_file(u[k]          , "log_alpha_f.asc", 	dim);
		write_single_file(u[k] + 1 * dim, "log_psi_f.asc",	dim);
		write_single_file(u[k] + 2 * dim, "phi_f.asc", 	dim);
		write_single_file(&w, "w_f.asc", 1);

		// Print final update.
		if (k > 0)
		{
			write_single_file(du[k - 1]          , "du0_f.asc", dim);
			write_single_file(du[k - 1] +     dim, "du1_f.asc", dim);
			write_single_file(du[k - 1] + 2 * dim, "du2_f.asc", dim);
		}

		// Print final RHS.
		write_single_file(f[k]          , "f0_f.asc", dim);
		write_single_file(f[k] +     dim, "f1_f.asc", dim);
		write_single_file(f[k] + 2 * dim, "f2_f.asc", dim);

		// Also print Newton parameters.
		switch (solverType)
		{
			case 1:
				write_single_file(norm_du,		"norm_du.asc",	 	k);
				write_single_file(norm_du_bar,	"norm_du_bar.asc",	k);
				break;
			case 2:
				write_single_file(norm_f,	"norm_f.asc",	 	k);
				break;
		}

		write_single_file(lambda,		"lambda.asc",	 k);
		write_single_file(Theta,		"Theta.asc",	 k);
		write_single_file(mu,		"mu.asc",	 k);
		write_single_file(lambda_prime,	"lambda_prime.asc",	 k);
		write_single_file(mu_prime,	"mu_prime.asc",		 k);

		// Print final iteration's RHS's norms.
		f_norms[0] = norm2(f[k]          );
		f_norms[1] = norm2(f[k] +     dim);
		f_norms[2] = norm2(f[k] + 2 * dim);
		printf("***                                                \n");
		printf("***        FINAL ITERATION:                        \n");
		printf("***           || f0 ||   = %-12.10E           \n", f_norms[0]);
		printf("***           || f1 ||   = %-12.10E           \n", f_norms[1]);
		printf("***           || f2 ||   = %-12.10E           \n", f_norms[2]);
		printf("***                                                \n");
		printf("******************************************************\n");

		// Also print omega.
		printf("******************************************************\n");
		printf("***                                                \n");
		printf("***           FINAL OMEGA:                         \n");
		printf("***            w          = %-12.10E            \n", w);
		printf("***                                                \n");
		printf("******************************************************\n");

		// ANALYSIS PHASE.
		ex_analysis(1, &M_ADM, &M_KOM, u[k], Dr_u, ghost, order, dim, dr);
		ex_phi_analysis(1, &phi_max, &hwl_res, u[k], ghost, order, dim);

		// Exit directory by going up one level (executable level).
		if ((io_code = chdir(work_dirname)) == -1)
		{
			printf("ERROR: Could not cd to output directory %s!\nError code = %d.\n", work_dirname, io_code);
			printf("errno : %s\n", strerror(errno));
			exit(-1);
		}

		// Rename directory to include w.
		snprintf(final_dirname, MAX_STR_LEN, "l=0,phi=%.5E,w=%.5E,dr=%.5E,N=%04lld,order=%lld", phi_max, w, dr, NrInterior, order);
		printf("Changing to directory %s.\n", final_dirname);
		if ((io_code = rename(initial_dirname, final_dirname)) != 0)
		{
			printf("ERROR: Could not rename directory %s into %s!\nError code = %d.\n", initial_dirname, final_dirname, io_code);
			printf("errno : %s\n", strerror(errno));
			exit(-1);
		}

		// Sweep continuation if sanity checks first.
		if (errCode == 0)
		{
			// Check if sweep should continue on this resolution.
			if (w <= w_min || w >= w_max)
			{
				printf("******************************************************\n");
				printf("***                                                \n");
				printf("***   Sweep cannot continue because w is out of range (%.5E, %.5E) !\n", w_min, w_max);
				printf("***                                                \n");
				printf("******************************************************\n");
				break;
			}
			else if (hwl_res < hwl_min)
			{
				printf("******************************************************\n");
				printf("***                                                \n");
				printf("***   Sweep cannot continue because N(HWL) < MIN(N(HWL)) = %lld !\n", hwl_min);
				printf("***   In other words, scalar field has not enough resolution. Try with more resolution or decrease hwl_min.\n");
				printf("***                                                \n");
				printf("******************************************************\n");
				break;
			}
			else if (hwl_res > hwl_max)
			{
				printf("******************************************************\n");
				printf("***                                                \n");
				printf("***   Sweep cannot continue because N(HWL) > MAX(N(HWL)) = %lld !\n", hwl_max);
				printf("***   In other words, scalar field is too scattered. Try with less resolution or increase hwl_max.\n");
				printf("***                                                \n");
				printf("******************************************************\n");
				break;
			}
			else
			{
				for (counter_i = 0; counter_i < GNUM; ++counter_i)
				{
					// Get peaks.
					peak_prev[counter_i] = u_seed[counter_i * dim + cblas_idamax(dim, u[k] + counter_i * dim, 1)];
					peak_next[counter_i] =   u[k][counter_i * dim + cblas_idamax(dim, u[k] + counter_i * dim, 1)];
				}
				next_scale[2] = peak_next[2] / peak_prev[2];
				next_scale[0] = 1.0 + next_scale[2] * (1.0 - peak_prev[0] / peak_next[0]);
				next_scale[1] = 1.0 + next_scale[2] * (1.0 - peak_prev[1] / peak_next[1]);

				for (counter_i = 0; counter_i < GNUM; ++counter_i)
					printf("**** Variable %lld peak = % -.5E, previous peak = % -.5E : predicted scale factor = %.5E\n", counter_i, peak_next[counter_i], peak_prev[counter_i], next_scale[counter_i]);
				
				// Omega prediction.
				peak_prev[GNUM] = omega_calc(u_seed[GNUM * dim], m);
				peak_next[GNUM] = omega_calc(  u[k][GNUM * dim], m);

				next_scale[GNUM] = 1.0 + next_scale[2] * (1.0 - peak_prev[GNUM] / peak_next[GNUM]);
				
				printf("**** scaled w = %.5E, w = %.5E, scale_u6 = %.5E\n", next_scale[GNUM] * w, w, next_scale[GNUM]);

				// Trasfer to initial data.
#ifdef NEXT_SCALE_JUMP
				#pragma omp parallel shared(u)
				{
					#pragma omp for schedule(dynamic, 1)
					for (counter_i = 0; counter_i < GNUM * dim + 1; ++counter_i)
					{
						u[0][counter_i] = -scale_next * u_seed[counter_i];
						u[0][counter_i] += (1.0 + scale_next) * u[k][counter_i];
						u_seed[counter_i] = u[k][counter_i];
					}
				}
#else
				#pragma omp parallel shared(u)
				{
					#pragma omp for schedule(dynamic, 1)
					for (counter_i = 0; counter_i < GNUM * dim + 1; ++counter_i)
					{
						u[0][counter_i] = u_seed[counter_i] = u[k][counter_i];
					}
				}
				// Scale variables.
				for (counter_i = 0; counter_i < GNUM; ++counter_i)
				{
					cblas_dscal(dim, next_scale[counter_i], u[0] + counter_i * dim, 1);
				}
				u[0][GNUM * dim] = inverse_omega_calc(next_scale[GNUM] * w, m);
#endif
				if (w_step != 0.0)
				{
					u[0][w_idx] = inverse_omega_calc(w + w_step, m);
				}
				// Set initial omega.
				w0 = omega_calc(u[0][w_idx], m);

				// Set analysis phase to 0 again.
				J.analysis_phase = 0;

				// Set new lambda0 to one since convergence has improved.
				lambda0 = 1.0;

				printf("******************************************************\n");
				printf("***                                                \n");
				printf("***   Setting initial data to last solution, scaling, and continuing...\n");
				printf("***                                                \n");
				printf("***                                                \n");
				printf("******************************************************\n");
			}
		}
		else
		{
			// Cannot continue sweep because errCode != 0.
			printf("******************************************************\n");
			printf("***                                                \n");
			printf("***   Sweep cannot continue because errCode = %lld !\n", errCode);
			printf("***                                                \n");
			printf("******************************************************\n");
			break;
		}
	}	

	// Print entire history of solutions and updates.

	// Clear memory.
	printf("***************************************************\n");
	printf("***                                                \n");
	printf("***              Deallocating memory...            \n");
	printf("***                                                \n");

	pardiso_stop();
	csr_deallocate(&J);

	for (i = 0; i < maxNewtonIter; i++)
	{
		SAFE_FREE(u[i]);
		SAFE_FREE(f[i]);
		SAFE_FREE(du[i]);
		SAFE_FREE(du_bar[i]);
	}

	SAFE_FREE(r);

	SAFE_FREE(u);
	SAFE_FREE(f);
	SAFE_FREE(du);
	SAFE_FREE(du_bar);

	SAFE_FREE(Dr_u);
	SAFE_FREE(Drr_u);

	SAFE_FREE(norm_f);
	SAFE_FREE(norm_du);
	SAFE_FREE(norm_du_bar);
	SAFE_FREE(lambda);
	SAFE_FREE(Theta);
	SAFE_FREE(mu);
	SAFE_FREE(lambda_prime);
	SAFE_FREE(mu_prime);
	
	SAFE_FREE(u_seed);

	// Clear libconfig configuration.
	config_destroy(&cfg);

	printf("***              Finished deallocation!            \n");
	printf("***                                                \n");
	printf("***************************************************\n");

	// Print final message.
	printf("***************************************************\n");
	printf("***                                                \n");
	printf("***           All done! Have a nice day!           \n");
	printf("***                                                \n");
	printf("***************************************************\n");
	printf("***************************************************\n");

	// All done.
	return 0;
}

// Centered print line.
//printf("***                        **                      \n");
