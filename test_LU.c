#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"
#include "f2c.h"

void print_permutation(integer n, int *p);
void print_matrix(integer n, doublereal *a);
void vector_copy(integer n, doublereal *u, doublereal *v);

int main(int argc, const char* argv[])
{
    FILE *fp;
    char *mode = "r";
    fp = fopen("donnees_LU.txt",mode);
    int n_mat;
    fscanf(fp,"%d",&n_mat);
    int i,j;
    integer n1, n2, n3;
    fscanf(fp,"%ld %ld %ld", &n1, &n2, &n3);
    doublereal* A;
    integer info;
    int* ipiv;

    printf("First Matrix\n");
    A = calloc(n1*n1,sizeof(doublereal));
    ipiv = calloc(n1,sizeof(int));
    for (i = 0; i < n1; i++) {
        fscanf(fp,"%lf %lf %lf", &(A[i]), &(A[i+1*n1]), &(A[i+2*n1]));
    }
    dgetrf_(&n1, &n1, A, &n1, ipiv, &info);
    print_permutation(n1, ipiv);
    print_matrix(n1, A);
    free(A);
    free(ipiv);

    printf("Second Matrix\n");
    A = calloc(n2*n2,sizeof(doublereal));
    ipiv = calloc(n2,sizeof(int));
    for (i = 0; i < n2; i++) {
        fscanf(fp,"%lf %lf %lf %lf", &(A[i]), &(A[i+1*n2]), &(A[i+2*n2]), &(A[i+3*n2]));
    }
    dgetrf_(&n2, &n2, A, &n2, ipiv, &info);
    print_permutation(n2, ipiv);
    print_matrix(n2, A);
    free(A);
    free(ipiv);

    printf("Third Matrix\n");
    A = calloc(n2*n2,sizeof(doublereal));
    ipiv = calloc(n2,sizeof(int));
    for (i = 0; i < n2; i++) {
        fscanf(fp,"%lf %lf %lf %lf", &(A[i]), &(A[i+1*n2]), &(A[i+2*n2]), &(A[i+3*n2]));
    }
    dgetrf_(&n2, &n2, A, &n2, ipiv, &info);
    print_permutation(n2, ipiv);
    print_matrix(n2, A);

    close(fp);
    return 0;
}

void print_permutation(integer n, int* p)
{
    int i,j;
    for (i = 0; i < n; i++) {
        p[i] = p[i]-1;
    }
    for (i = 0; i < n; i++) {
        p[p[i]] = i;
        for (j = 0; j < n; j++) {
            printf("%d ", j==p[i]);
        }
        printf("\n");
    }

    printf("\n");
}

void print_matrix(integer n, doublereal *a)
{
    integer i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", a[i+j*n]);
        }
        printf("\n");
    }
    printf("\n");
}

// copy vector u to vector v
void vector_copy(integer n, doublereal *u, doublereal *v)
{
    integer i;
    for (i = 0; i < n; i++) {
        v[i] = u[i];
    }
}
