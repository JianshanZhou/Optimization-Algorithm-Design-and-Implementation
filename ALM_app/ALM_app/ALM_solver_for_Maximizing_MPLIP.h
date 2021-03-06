/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ALM_solver_for_Maximizing_MPLIP.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 30-May-2019 11:43:08
 */

#ifndef ALM_SOLVER_FOR_MAXIMIZING_MPLIP_H
#define ALM_SOLVER_FOR_MAXIMIZING_MPLIP_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "ALM_solver_for_Maximizing_MPLIP_types.h"

/* Function Declarations */
#ifdef __cplusplus

extern "C" {

#endif

  extern void ALM_solver_for_Maximizing_MPLIP(double amp, const double phase[3],
    double relFreqConst, double numSources, const double xcoords[3], const
    double ycoords[3], const double zcoords[3], double sigma, double epsilon,
    double alpha, double beta, double maxIter, const double lb[2], const double
    ub[2], double x[2], double *fval);

#ifdef __cplusplus

}
#endif
#endif

/*
 * File trailer for ALM_solver_for_Maximizing_MPLIP.h
 *
 * [EOF]
 */
