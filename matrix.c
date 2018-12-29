#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void print_vector(int n, double *x)
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%le ", x[i]);
    }
    printf("\n\n");
}

void print_matrix(int m, int n, double *a)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%le ", a[i+j*m]);
        }
        printf("\n");
    }
    printf("\n");
}


void read_vector(FILE* fp, int n, double *v) {
    char line[1024], *p, *e;
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
    free(v);
}

// read a symetric matrix in column form
void read_matrix_sc(FILE* fp, int n, double *a_vec) {
    int i,j;
    double val;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fscanf(fp,"%lf\n", &val);
            if (j <= i) {
                a_vec[j*n + i - (j*(j+1))/2] = val;
            }
        }
    }
}

void write_vector_c(FILE *fp, int n, double *v) {
    int i;
    fprintf(fp, "%d\n",n);
    for (i = 0; i < n; i++) {
       fprintf(fp, "%le\n",v[i]); 
    }
}

void write_vector(FILE *fp, int n, double *v) {
    int i;
    for (i = 0; i < n; i++) {
        fprintf(fp, "%le ", v[i]);
    }
    fprintf(fp, "\n");
}

void copy_vector(int n, double *u, double *v) {
    int i;
    for (i = 0; i < n; i++) {
        v[i] = u[i];
    }
}

void gaxpy_s(int n, double *a_vec, double *x, double *y) {
    int i,j;

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            y[i] += a_vec[i*n + j - (i*(i+1))/2]*x[j];
        }
        for (i = j; i < n; i++) {
            y[i] += a_vec[j*n + i - (j*(j+1))/2]*x[j];
        }
    }
}

void saxpy(int n, double a, double *x, double *y) {
    int i;
    for (i = 0; i < n; i++) {
        y[i] += a*x[i];
    }
}

void vector_nullify(int n, double *v) {
    int i;
    for (i = 0; i < n; i++) {
        v[i] = 0.0;
    }
}

// computes the dot product using Kahan sum algorithm
double vector_dot_product(int n, double *u, double *v) {
    double sum = 0.0;
    double c = 0.0;
    double y,t;
    int i;

    for (i = 0; i < n; i++) {
        y = u[i]*v[i] - c;
        t = sum + y;
        c = (t-sum)-y;
        sum = t;
    }
    return sum;
}

double vector_norm(int n, double *v) {
    return sqrt(vector_dot_product(n,v,v));
}
