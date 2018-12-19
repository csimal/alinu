#include <stdio.h>
#include <stdlib.h>

#ifndef MATRIX_H
#define MATRIX_H

void print_vector(int n, double *v);

void print_matrix(int m, int n, double *a);

void copy_vector(int n, double *u, double *v);

void read_vector(FILE* fp, int n, double *v);

void read_matrix(FILE* fp, int m, int n, double *a);


#endif
