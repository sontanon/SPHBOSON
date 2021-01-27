#include "tools.h"
#include "param.h"

// Macros for parameter ranges.
#define MAX_DR 1.0
#define MIN_DR 0.001

#define MAX_NRINTERIOR 9999999LL
#define MIN_NRINTERIOR 32LL

#define MAX_M 1.0E+3
#define MIN_M 1.0E-3

#define MAX_PHI0 1.0E+5
#define MIN_PHI0 1.0E-5

#define MAX_W0 1.0
#define MIN_W0 0.0

#define MAX_MAXITER 100000LL
#define MIN_MAXITER 0LL

#define MIN_WEIGHT 1.0E-16

#define MAX_EPS 1.0E-1
#define MIN_EPS 1.0E-16

void parser(const char *fname)
{
	// Initialize cfg.
	config_init(&cfg);

	// Read the file. If there is an error, report and exit.
	if (!config_read_file(&cfg, fname))
	{
		fprintf(stderr, "PARSER: CRITICAL ERROR IN FILE!\n");
		fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
		config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		exit(-1);
	}

	// Parse arguments doing sanity checks.

	// GRID.
	// dr.
	if (config_lookup_float(&cfg, "dr", &dr) == CONFIG_TRUE)
	{
		if (MAX_DR < dr || dr < MIN_DR)
		{
			fprintf(stderr, "PARSER: ERROR! dr = %3.5E is not in range [%3.5E, %3.5E]\n", dr, MIN_DR, MAX_DR);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"dr\" from parameter file. Setting to default value, dr = %3.5E\n", dr);
	}
	// NrInterior.
	if (config_lookup_int64(&cfg, "NrInterior", &NrInterior) == CONFIG_TRUE)
	{
		if (MAX_NRINTERIOR < NrInterior || NrInterior < MIN_NRINTERIOR)
		{
			fprintf(stderr, "PARSER: ERROR! NrInterior = %lld is not in range [%lld, %lld]\n", NrInterior, MIN_NRINTERIOR, MAX_NRINTERIOR);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"NrInterior\" from parameter file. Setting to default value, NrInterior = %lld\n", NrInterior);
	}
	// order.
	if (config_lookup_int64(&cfg, "order", &order) == CONFIG_TRUE)
	{
		if (order != 2 && order != 4)
		{
			fprintf(stderr, "PARSER: ERROR! order = %lld is not supported. Only 2 or 4 are supported finite difference orders.\n", order);
			fprintf(stderr, "        Please input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"order\" from parameter file. Setting to default value, order = %lld\n", order);
	}
	// Set ghost zones.
	ghost = order / 2;

	// Do not forget to change dimension parameter!
	dim = NrInterior + 2 * ghost;
	w_idx = GNUM * dim;

	// SCALAR FIELD PARAMETERS.
	// m.
	if (config_lookup_float(&cfg, "m", &m) == CONFIG_TRUE)
	{
		if (MAX_M < m || m < MIN_M)
		{
			fprintf(stderr, "PARSER: ERROR! m = %3.5E is not in range [%3.5E, %3.5E]\n", m, MIN_M, MAX_M);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"m\" value from parameter file. Setting to default value, m = %3.5E\n", m);
	}
	// phi0.
	if (config_lookup_float(&cfg, "phi0", &phi0) == CONFIG_TRUE)
	{
		if (MAX_PHI0 < phi0 || phi0 < MIN_PHI0)
		{
			fprintf(stderr, "PARSER: ERROR! phi0 = %3.5E is not in range [%3.5E, %3.5E]\n", phi0, MIN_PHI0, MAX_PHI0);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"phi0\" value from parameter file. Setting to default value, phi0 = %3.5E\n", phi0);
	}
	// sigma.
	config_lookup_float(&cfg, "sigma", &sigma);
	
	// w0.
	if (config_lookup_float(&cfg, "w0", &w0) == CONFIG_TRUE)
	{
		if (MAX_W0 < w0 / m || w0 / m < MIN_W0)
		{
			fprintf(stderr, "PARSER: ERROR! (w0 / m) = (%3.5E / m) is not in range (%3.5E, %3.5E)\n", w0, MIN_W0, MAX_W0);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"w0\" value from parameter file. Setting to default value, w0 = %3.5E\n", w0);
	}	
	// fixedPhi.
	if (config_lookup_int64(&cfg, "fixedPhi", &fixedPhi) == CONFIG_TRUE)
	{
		if (NrInterior < fixedPhi || fixedPhi < 0)
		{
			fprintf(stderr, "PARSER: ERROR! fixedPhi = %lld is not in range [%lld, %lld]\n", fixedPhi, 0LL, NrInterior);
			fprintf(stderr, "        Please input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"fixedPhi\" from parameter file. Setting to default value, fixedPhi = %lld\n", fixedPhi);
	}

	// INITIAL DATA.
	// readInitialData.
	if (config_lookup_int64(&cfg, "readInitialData", &readInitialData) == CONFIG_TRUE)
	{
		if (readInitialData != 0 && readInitialData != 1 && readInitialData != 2 && readInitialData != 3)
		{
			fprintf(stderr, "PARSER: ERROR! readInitialData = %lld is not supported. Only 0, 1, 2, or 3 as boolean values for indication of whether to read initial data specified by user.\n", readInitialData);
			fprintf(stderr, "        Please input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"readInitialData\" from parameter file. Setting to default value, readInitialData = %lld\n", readInitialData);
	}
	// Read initial data parameters.
	switch (readInitialData)
	{
		// Interpolation from different size and/or resolution grid.
		case 3:
			// Read filenames.
			if (config_lookup_string(&cfg, "log_alpha_i", &log_alpha_i) == CONFIG_FALSE)
			{
				fprintf(stderr, "PARSER: ERROR! readInitialData = 3 requires values for all initial files. Did not find \"log_alpha_i\".\n");
				exit(-1);
			}
			if (config_lookup_string(&cfg, "log_psi_i", &log_psi_i) == CONFIG_FALSE)
			{
				fprintf(stderr, "PARSER: ERROR! readInitialData = 3 requires values for all initial files. Did not find \"log_psi_i\".\n");
				exit(-1);
			}
			if (config_lookup_string(&cfg, "phi_i", &phi_i) == CONFIG_FALSE)
			{
				fprintf(stderr, "PARSER: ERROR! readInitialData = 3 requires values for all initial files. Did not find \"phi_i\".\n");
				exit(-1);
			}

			// Grid parameters.
			if (config_lookup_int64(&cfg, "dimInitial", &dimInitial) == CONFIG_TRUE)
			{
				if (MAX_NRINTERIOR < dimInitial || dimInitial < MIN_NRINTERIOR)
				{
					fprintf(stderr, "PARSER: ERROR! dimInitial = %lld is not in range [%lld, %lld]\n", dimInitial, MIN_NRINTERIOR, MAX_NRINTERIOR);
					fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
					exit(-1);
				}
			}
			else
			{
				fprintf(stderr, "PARSER: ERROR! readInitialData = 3 requires value for dimInitial.\n");
				exit(-1);
			}
			if (config_lookup_int64(&cfg, "order_i", &order_i) == CONFIG_TRUE)
			{
				if (order_i != 2 && order_i != 4)
				{
					fprintf(stderr, "PARSER: ERROR! order_i = %lld is not supported. Only 2 or 4 are supported finite difference orders.\n", order);
					fprintf(stderr, "        Please input proper value in parameter file.\n");
					exit(-1);
				}
			}
			else
			{
				fprintf(stderr, "PARSER: ERROR! readInitialData = 3 requires value for order_i.\n");
				exit(-1);
			}
			if (config_lookup_int64(&cfg, "ghost_i", &ghost_i) == CONFIG_TRUE)
			{
				if (ghost_i != 1 && ghost_i != 2)
				{
					fprintf(stderr, "PARSER: ERROR! ghost_i = %lld is not supported. Only 1 or 2 are supported.\n", order);
					fprintf(stderr, "        Please input proper value in parameter file.\n");
					exit(-1);
				}
			}
			else
			{
				fprintf(stderr, "PARSER: ERROR! readInitialData = 3 requires value for ghost_i.\n");
				exit(-1);
			}			
			if (config_lookup_float(&cfg, "dr_i", &dr_i) == CONFIG_TRUE)
			{
				if (MAX_DR < dr_i || dr_i < MIN_DR)
				{
					fprintf(stderr, "PARSER: ERROR! dr_i = %3.5E is not in range [%3.5E, %3.5E]\n", dr_i, MIN_DR, MAX_DR);
					fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
					exit(-1);
				}
			}
			else
			{
				fprintf(stderr, "PARSER: ERROR! readInitialData = 3 requires value for dr_i.\n");
				exit(-1);
			}
			config_lookup_string(&cfg, "w_i", &w_i);

			break;
	
		// Default case for 1 or 2.
		case 2:
		case 1:
			config_lookup_string(&cfg, "log_alpha_i", &log_alpha_i);
			config_lookup_string(&cfg, "log_psi_i", &log_psi_i);
			config_lookup_string(&cfg, "phi_i", &phi_i);
			
			// Initial Data extensions.
			if (readInitialData == 2)
			{
				// NrTotalInitial.
				config_lookup_int64(&cfg, "dimInitial", &dimInitial);
			}
			else
			{
				dimInitial = dim;
			}
			break;
	}

	// Scale initial data.
	config_lookup_float(&cfg, "scale_u0", &scale_u0);
	config_lookup_float(&cfg, "scale_u1", &scale_u1);
	config_lookup_float(&cfg, "scale_u2", &scale_u2);
	config_lookup_float(&cfg, "scale_u3", &scale_u3);

	// Generate via analytic guess.
	if (!readInitialData)
	{
		// psi0.
		if (config_lookup_float(&cfg, "phi0", &phi0) == CONFIG_TRUE)
		{
			if (MAX_PHI0 < phi0 || phi0 < MIN_PHI0)
			{
				fprintf(stderr, "PARSER: ERROR! phi0 = %3.5E is not in range [%3.5E, %3.5E]\n", phi0, MIN_PHI0, MAX_PHI0);
				fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
				exit(-1);
			}
		}
		else
		{
			fprintf(stderr, "PARSER: WARNING! Could not properly read \"phi0\" value from parameter file. Setting to default value, phi0 = %3.5E\n", phi0);
		}
	}

	// Initial frequency.
	if (!w_i)
	{
		// w0.
		if (config_lookup_float(&cfg, "w0", &w0) == CONFIG_TRUE)
		{
			if (MAX_W0 < w0 / m || w0 / m < MIN_W0)
			{
				fprintf(stderr, "PARSER: ERROR! (w0 / m) = (%3.5E / m) is not in range (%3.5E, %3.5E)\n", w0, MIN_W0, MAX_W0);
				fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
				exit(-1);
			}
		}
		else
		{
			fprintf(stderr, "PARSER: WARNING! Could not properly read \"w0\" value from parameter file. Setting to default value, w0 = %3.5E\n", w0);
		}	
	}

	// SOLVER PARAMETERS.
	// solverType.
	if (config_lookup_int64(&cfg, "solverType", &solverType) == CONFIG_TRUE)
	{
		if (solverType != 1 && solverType != 2 && solverType != 3)
		{
			fprintf(stderr, "PARSER: ERROR! solverType = %lld is not supported. Only 1 or 2 or 3 are supported for indication of whether to use error or residual based solver.\n", solverType);
			fprintf(stderr, "        Please input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"solverType\" from parameter file. Setting to default value, solverType = %lld\n", solverType);
	}
	// localSolver.
	if (config_lookup_int64(&cfg, "localSolver", &localSolver) == CONFIG_TRUE)
	{
		if (localSolver != 0 && localSolver != 1)
		{
			fprintf(stderr, "PARSER: ERROR! localSolver = %lld is not supported. Only 0 or 1 boolean are supported for indication of whether to use local solver inside global solver.\n", localSolver);
			fprintf(stderr, "        Please input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"localSolver\" from parameter file. Setting to default value, localSolver = %lld\n", localSolver);
	}
	// epsilon.
	if (config_lookup_float(&cfg, "epsilon", &epsilon) == CONFIG_TRUE)
	{
		if (MAX_EPS < epsilon || epsilon < MIN_EPS)
		{
			fprintf(stderr, "PARSER: ERROR! exit tolerance epsilon = %3.5E is not in range [%3.5E, %.35E]\n", epsilon, MIN_EPS, MAX_EPS);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"epsilon\" from parameter file. Setting to default value, epsilon = %3.5E\n", epsilon);
	}
	// maxNewtonIter.
	if (config_lookup_int64(&cfg, "maxNewtonIter", &maxNewtonIter) == CONFIG_TRUE)
	{
		if (MAX_MAXITER < maxNewtonIter || maxNewtonIter < MIN_MAXITER)
		{
			fprintf(stderr, "PARSER: ERROR! maxNewtonIter = %lld is not in range [%lld, %lld]\n", maxNewtonIter, MIN_MAXITER, MAX_MAXITER);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"maxNewtonIter\" from parameter file. Setting to default value, maxIter = %lld\n", maxNewtonIter);
	}
	// lambda0.
	if (config_lookup_float(&cfg, "lambda0", &lambda0) == CONFIG_TRUE)
	{
		if (1.0 < lambda0 || lambda0 <= MIN_WEIGHT)
		{
			fprintf(stderr, "PARSER: ERROR! initial damping factor lambda0 = %3.5E is not in range (%3.5E, 1.0]\n", lambda0, MIN_WEIGHT);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"lambda0\" from parameter file. Setting to default value, lambda0 = %3.5E\n", lambda0);
	}
	// lambdaMin.
	if (config_lookup_float(&cfg, "lambdaMin", &lambdaMin) == CONFIG_TRUE)
	{
		if (lambda0 <= lambdaMin || lambdaMin < MIN_WEIGHT)
		{
			fprintf(stderr, "PARSER: ERROR! minimum damping factor lambdaMin = %3.5E is not in range [%3.5E, lambda0)\n", lambdaMin, MIN_WEIGHT);
			fprintf(stderr, "        Please edit range in \"parser.c\" source file or input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"lambdaMin\" from parameter file. Setting to default value, lambdaMin = %3.5E\n", lambdaMin);
	}
	// useLowRank.
	if (config_lookup_int64(&cfg, "useLowRank", &useLowRank) == CONFIG_TRUE)
	{
		if (useLowRank != 0 && useLowRank != 1)
		{
			fprintf(stderr, "PARSER: ERROR! useLowRank = %lld is not supported. Only 0 or 1 boolean are supported for indication of whether to use Low Rank Update.\n", useLowRank);
			fprintf(stderr, "        Please input proper value in parameter file.\n");
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"useLowRank\" from parameter file. Setting to default value, useLowRank = %lld\n", useLowRank);
	}

	// INITIAL GUESS CHECK.
	if (config_lookup_int64(&cfg, "max_initial_guess_checks", &max_initial_guess_checks) == CONFIG_TRUE)
	{
		if (max_initial_guess_checks < 0 || max_initial_guess_checks > 10)
		{
			fprintf(stderr, "PARSER: ERROR! max_initial_guess_checks = %lld is out of bounds.\n", max_initial_guess_checks);
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"max_initial_guess_checks\" from parameter file. Setting to default value, max_initial_guess_checks = %lld\n", max_initial_guess_checks);
	}
	if (config_lookup_float(&cfg, "norm_f0_target", &norm_f0_target) == CONFIG_TRUE)
	{
		if (norm_f0_target < epsilon || norm_f0_target > 1.0)
		{
			fprintf(stderr, "PARSER: ERROR! norm_f0_target = %3.5E out of bounds!\n", norm_f0_target);
			exit(-1);
		}
	}
	else
	{
		fprintf(stderr, "PARSER: WARNING! Could not properly read \"norm_f0_target\" from parameter file. Setting to default value, norm_f0_target = %3.5E\n", norm_f0_target);
	}
	
	// SWEEP CONTROL.
	config_lookup_int64(&cfg, "hwl_min", &hwl_min);
	config_lookup_int64(&cfg, "hwl_max", &hwl_max);
	config_lookup_float(&cfg, "w_max", &w_max);
	config_lookup_float(&cfg, "w_min", &w_min);
	config_lookup_float(&cfg, "w_step", &w_step);

	// NEXT SCALE ADVANCEMENT.
	config_lookup_float(&cfg, "scale_next", &scale_next);

	// OUTPUT
	// work_dirname.
	getcwd(work_dirname, MAX_STR_LEN);

	// Set initial directory name.
	snprintf(initial_dirname, MAX_STR_LEN, "l=0,phi=X.XXXXXE+00,w=X.XXXXXE-01,dr=%.5E,N=%04lld,order=%lld", dr, NrInterior, order);

	// All done.
	return;
}
