/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ALM_solver_for_Maximizing_MPLIP_emxutil.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 30-May-2019 11:43:08
 */

#ifndef ALM_SOLVER_FOR_MAXIMIZING_MPLIP_EMXUTIL_H
#define ALM_SOLVER_FOR_MAXIMIZING_MPLIP_EMXUTIL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "ALM_solver_for_Maximizing_MPLIP_types.h"

/* Function Declarations */
#ifdef __cplusplus

extern "C" {

#endif

  extern void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel);
  extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
  extern void emxFree_char_T(emxArray_char_T **pEmxArray);
  extern void emxFree_real_T(emxArray_real_T **pEmxArray);
  extern void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);
  extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#ifdef __cplusplus

}
#endif
#endif

/*
 * File trailer for ALM_solver_for_Maximizing_MPLIP_emxutil.h
 *
 * [EOF]
 */
