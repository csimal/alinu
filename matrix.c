#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>

void print_vector(int n, double *x)
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%e ", x[i]);
    }
    printf("\n\n");
}

void print_matrix(int m, int n, double *a)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%e ", a[i+j*m]);
        }
        printf("\n");
    }
    printf("\n");
}

// copy vector u to vector v
void copy_vector(int n, double *u, double *v)
{
    int i;
    for (i = 0; i < n; i++) {
        v[i] = u[i];
    }
}

void read_vector(FILE* fp, int n, double *v) {
    char line[1024], *p, *e;
    double val;
    int i;

    fgets(line, sizeof(line), fp);
    for (i = 0, p = line; i < n; i++, p = e) {
        v[i] = strtod(p, &e);
    }
}

// read matrix from file and put it in column major order in a
void read_matrix(FILE* fp, int m, int n, double *a) {
    double *v = calloc(n, sizeof(double));
    int i, j;

    for (i = 0; i < m; i++) {
        read_vector(fp, n, v);
        for (j = 0; j < n; j++) {
            a[i+j*m] = v[j];
        }
    }
}

void read_matrix_2(FILE* fp, int m, int n, double* a) {
    int i,j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            fscanf(fp,"%lf\n", &(a[i+j*m]));
        }
    }
}
