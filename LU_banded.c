#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "f2c.h"
#include "matrix.h"

#define MIN(X, Y) (((X)<(Y))? (X):(Y))
#define MAX(X, Y) (((X)>(Y))? (X):(Y))

void print_permutation(int n, int *p);

int main(int argc, const char* argv[])
{
    FILE *fp;
    int n,i,j;
    doublereal *A;
    int *ipiv;
    int info;

    fp = fopen("Matrice_bande.txt","r");
    fscanf(fp,"%d\n",&n);
    printf("%d\n",n);
    A = calloc(n*n,sizeof(doublereal));
    read_matrix(fp,n,n,A);
    printf("A:\n");
    print_matrix(n,n,A);
    ipiv = calloc(n,sizeof(int));
    integer Kl = 2, Ku = 3, ldab = 2*Kl+Ku+1;
    doublereal *AB = calloc(ldab*n,sizeof(doublereal));

    for (j = 0; j < n; j++) {
        for (i = MAX(0,j-Ku); i < MIN(n,j+Kl+1); i++) {
            AB[Kl+Ku+i-j+j*ldab] = A[i+j*n];
        }
    }
    fclose(fp);

    dgbtrf_(&n,&n,&Kl,&Ku,AB,&ldab,ipiv,&info);

    printf("info: %d\n",info);

    doublereal *L,*U;
    U = calloc(n*n,sizeof(doublereal));
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            if (i>=j-Kl-Ku)
                U[i+j*n] = AB[i-j+Kl+Ku+j*ldab];
        }
    }
    L = calloc(n*n,sizeof(doublereal));
    for (j = 0; j < n; j++) {
        L[j+j*n] = 1.0;
        int mi = MIN(n,j+1+Kl);
        for (i = j+1; i < mi; i++) {
            L[i+j*n] = AB[i-j+Ku+Kl+j*ldab];
        }
    }
    for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            doublereal perm = L[i+j*n];
            L[i+j*n] = L[ipiv[i]-1+j*n];
            L[ipiv[i]-1+j*n] = perm;
        }
    }
    printf("Raw result:\n");
    print_matrix(ldab,n,AB);
    printf("U:\n");
    print_matrix(n,n,U);
    printf("L:\n");
    print_matrix(n,n,L);
    printf("P:\n");
    print_permutation(n,ipiv);

    free(L);
    free(U);
    free(ipiv);
    free(A);

    return 0;
}

void print_permutation(int n, int* p)
{
    int i,j,perm;
    int *P = calloc(n*n,sizeof(int));
    for (i = 0; i < n; i++) {
        P[i+i*n] = 1;
        p[i] = p[i]-1;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            perm = P[i+j*n];
            P[i+j*n] = P[p[i]+j*n];
            P[p[i]+j*n] = perm;
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%d ",P[i+j*n]);
        }
        printf("\n");
    }
    free(P);
}
