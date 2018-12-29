#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "matrix.h"

#ifndef CG_H
#define CG_H

int solve_cg(FILE *fp, int n, double *a_vec, double *b, double *x, double tol);

int solve_cg_jacobi(FILE *fp, int n, double *a_vec, double *b, double *x, double tol);

int solve_cg_ssor(FILE *fp, int n, double *a_vec, double *b, double *x, double w, double tol);

int solve_cg_cholesky(FILE *fp, int n, double *a_vec, double *b, double *x, double *k_vec, double tol);

int solve_cg_spectral(FILE *fp, int n, double *a_vec, double *b, double *x, int m, double tol);

double* eig(int n, double *a_vec, int m, char jobz, double *v);

#endif
