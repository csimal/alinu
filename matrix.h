#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"

#ifndef MATRIX_H
#define MATRIX_H

void print_vector(int n, doublereal *v);

void print_matrix(int m, int n, doublereal *a);

void copy_vector(int n, doublereal *u, doublereal *v);

void read_vector(FILE* fp, int n, doublereal *v);

void read_matrix(FILE* fp, int m, int n, doublereal *a);


#endif
