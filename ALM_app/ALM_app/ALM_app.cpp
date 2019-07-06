// ALM_app.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include <time.h>
#include <math.h>

/* Include Files */
#include "rt_nonfinite.h"
#include "ALM_solver_for_Maximizing_MPLIP.h"
#include "ALM_solver_for_Maximizing_MPLIP_terminate.h"
#include "ALM_solver_for_Maximizing_MPLIP_initialize.h"

#define PI 3.14159265358979323846
#define CLOCKS_PER_SEC ((clock_t)1000)

using namespace std;

/* Function Declarations */
static void argInit_2x1_real_T(double result[2]);
static void argInit_3x1_real_T(double result[3]);
static double argInit_real_T(void);
static void main_ALM_solver_for_Maximizing_MPLIP(void);

/* Function Definitions */

/*
 * Arguments    : double result[2]
 * Return Type  : void
 */
static void argInit_2x1_real_T(double result[2])
{
	double result_tmp;

	/* Loop over the array to initialize each element. */
	/* Set the value of the array element.
	   Change this value to the value that the application requires. */
	result_tmp = argInit_real_T();
	result[0] = result_tmp;

	/* Set the value of the array element.
	   Change this value to the value that the application requires. */
	result[1] = result_tmp;
}

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
	double result_tmp;

	/* Loop over the array to initialize each element. */
	/* Set the value of the array element.
	   Change this value to the value that the application requires. */
	result_tmp = argInit_real_T();
	result[0] = result_tmp;

	/* Set the value of the array element.
	   Change this value to the value that the application requires. */
	result[1] = result_tmp;

	/* Set the value of the array element.
	   Change this value to the value that the application requires. */
	result[2] = argInit_real_T();
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
	return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_ALM_solver_for_Maximizing_MPLIP(void)
{
	double amp;
	double phase_tmp[3];
	double relFreqConst;
	double numSources;
	double xcoords[3];
	double ycoords[3];
	double zcoords[3];
	double sigma;
	double epsilon;
	double alpha;
	double beta;
	double maxIter;
	double lb[2];
	double ub[2];
	double x[2];
	double fval;

	/* Initialize function 'ALM_solver_for_Maximizing_MPLIP' input arguments. */
	amp = argInit_real_T();

	/* Initialize function input argument 'phase'. */
	argInit_3x1_real_T(phase_tmp);
	relFreqConst = argInit_real_T();
	numSources = argInit_real_T();

	/* Initialize function input argument 'xcoords'. */
	argInit_3x1_real_T(xcoords);

	/* Initialize function input argument 'ycoords'. */
	argInit_3x1_real_T(ycoords);

	/* Initialize function input argument 'zcoords'. */
	argInit_3x1_real_T(zcoords);
	sigma = argInit_real_T();
	epsilon = argInit_real_T();
	alpha = argInit_real_T();
	beta = argInit_real_T();
	maxIter = argInit_real_T();

	/* Initialize function input argument 'lb'. */
	argInit_2x1_real_T(lb);
	argInit_2x1_real_T(ub);

	/* Assign the initial values to the input arguments. */
	amp = 2.2;
	relFreqConst = (2 * PI) * 2.5;
	numSources = 3;
	double height = 3;
	xcoords[0] = 2.4112; xcoords[1] = 0.2064; xcoords[2] = 1.6787;
	ycoords[0] = 0.3957; ycoords[1] = 0.3927; ycoords[2] = 0.9877;
	zcoords[0] = height; zcoords[1] = height; zcoords[2] = height;
	lb[0] = -0.5; lb[1] = -2;
	ub[0] = 3.5; ub[1] = 3;

	sigma = 3;
	epsilon = 1e-6;
	alpha = 1.3;
	beta = 0.5;
	maxIter = 1000;

	long num_of_runs = 1000;
	double time_step = 0.05;
	clock_t start_T, end_T;

	/* Initialize function input argument 'ub'. */
	/* Call the entry-point 'ALM_solver_for_Maximizing_MPLIP'. */
	start_T = clock();
	for (long i = 0; i < num_of_runs; ++i)
	{
		ALM_solver_for_Maximizing_MPLIP(amp, phase_tmp, relFreqConst, numSources,
			xcoords, ycoords, zcoords, sigma, epsilon, alpha, beta, maxIter, lb,
			ub, x, &fval);
	}
	end_T = clock();

	/* Show the optimization result. */
	std::cout << "The optimal point is [" << x[0] << "," << x[1] << "]" << std::endl;
	std::cout << "The maximum pattern is " << -fval << std::endl;
	cout << "Avg. algorithm running time: " << ((double)(end_T - start_T) / CLOCKS_PER_SEC) / num_of_runs << "seconds." << endl;
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
	(void)argc;
	(void)argv;

	/* Initialize the application.
	   You do not need to do this more than one time. */
	ALM_solver_for_Maximizing_MPLIP_initialize();

	/* Invoke the entry-point functions.
	   You can call entry-point functions multiple times. */
	main_ALM_solver_for_Maximizing_MPLIP();

	/* Terminate the application.
	   You do not need to do this more than one time. */
	ALM_solver_for_Maximizing_MPLIP_terminate();
	return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */