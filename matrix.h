#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MATRIX_H
#define MATRIX_H

void print_vector(int n, double *v);

void print_matrix(int m, int n, double *a);

void read_vector(FILE* fp, int n, double *v);

void read_matrix(FILE* fp, int m, int n, double *a);

void read_matrix_sc(FILE* fp, int n, double *a_vec);

void write_vector_c(FILE *fp, int n, double *v);

void write_vector(FILE *fp, int n, double *v);

void copy_vector(int n, double *u, double *v);

void gaxpy_s(int n, double *a_vec, double *x, double *y);

void saxpy(int n, double a, double *x, double *y);

void vector_nullify(int n, double *v);

double vector_dot_product(int n, double *u, double *v);

double vector_norm(int n, double *v);

double* packed_to_full(int n, double *a_vec);

double* full_to_packed(int n, double *a_full);

#endif
