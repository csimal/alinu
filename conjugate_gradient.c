#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define TOL 0.001

double vector_norm(int n, double *v);
double* solve_cg(int n, double *a, double *b);
double* solve_pcg(int n, double *a, double *b, double *m);

int main(int argc, char const *argv[]) {
    FILE* fp;
    double *a, *b;
    int i,n;
    char filename[16];
    int lcm = LAPACK_COL_MAJOR;


    fp = fopen("donneeRHS.dat","r");
    fscanf(fp,"%d\n",&n);
    b = calloc(n,sizeof(double));
    a = calloc(n*n,sizeof(double));
    read_matrix_2(fp,1,n,b);
    fclose(fp);

    for (i = 1; i < 4; i++) {
        sprintf(filename,"donneeSyst%d.dat",i);
        fp = fopen(filename,"r");
        fscanf(fp,"%d\n",&n);
        read_matrix_2(fp,n,n,a);
    }



    return 0;
}

/*
* returns the euclidian norm of vector v
*/
double vector_norm(int n, double *v) {
    double sum = 0.0;
    int i;
    for (i = 0; i < n; i++)
        sum += v[i]*v[i];

    return sqrt(sum);
}

double* solve_cg(int n, double *a, double *b) {
    double *x = calloc(n,sizeof(double)); //NB calloc initialises the allocated memory to 0
    double *r = calloc(n,sizeof(double));
    double *p = calloc(n,sizeof(double));
    int i,j,k;

    for (i = 0; i < n; i++) {
        r[i] = -b[i];
        for (j = 0; j < n; j++) {
            r[i] += a[i+j*m]*x[i];
        }
        p[i] = r[i];
    }
    while (vector_norm(n,r) > TOL) {

    }
}

double* solve_pcg(int n, double *a, double *b, double *m) {
    return calloc(1,sizeof(double));
}
