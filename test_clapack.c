#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"
#include "f2c.h"
//#include <clapack.h>

void print_vector(integer n, doublereal *x);
void print_matrix(integer n, doublereal *a);
void vector_copy(integer n, doublereal *u, doublereal *v);

int main(int argc, const char* argv[])
{
    FILE *fp;
    char *mode = "r";
    fp = fopen("donnees_Lxb.txt",mode);

    integer i, j;
    integer n, nhrs, info, lwork;
    nhrs = 1;
    fscanf(fp,"%ld",&n);

    integer *ipvt = calloc(n,sizeof(int));
    doublereal *A_0, *A;
    doublereal *x, *b, *s;
    doublereal *work;
    lwork = 5*n;
    A = calloc(n*n, sizeof(double));
    A_0 = calloc(n*n, sizeof(double));
    b = calloc(n, sizeof(double));
    x = calloc(n,sizeof(double));
    s = calloc(n,sizeof(double));

    printf("alloc done\n");

    for (i = 0; i < n; i++) {
        fscanf(fp, "%le %le %le %le", &(A[i+0*n]), &(A[i+1*n]), &(A[i+2*n]), &(A[i+3*n]) );
    }
    vector_copy(n*n,A,A_0);
    printf("A initialized\n");
    for (i = 0; i < n; i++) {
        fscanf(fp, "%le", &(b[i]));
    }
    fclose(fp);

    vector_copy(n,b,x);

    printf("init done\n");

    print_matrix(n,A);
    print_vector(n,x);

    printf("solving system:\n");

    printf("solving with dtrtrs:\n");
    dtrtrs_("L","N","N", &n, &nhrs, A, &n, x, &n, &info);
    print_vector(n,x);
    printf("solving with dgesv:\n");
    vector_copy(n,b,x);
    dgesv_(&n, &nhrs, A, &n, ipvt, x, &n, &info);
    print_vector(n,x);

    vector_copy(n*n,A_0,A);

    printf("Condition numbers\n");
    work = calloc(5*n, sizeof(double));
    dgesvd_("N","N", &n, &n, A, &n, s, NULL, &n, NULL, &n, work, &lwork, &info);
    printf("%lf\n", s[0]/s[n-1]); //compute 2-norm condition number
    print_vector(n,s);

    vector_copy(n*n,A_0,A);
    doublereal rcond;
    free(work);
    work = calloc(3*n,sizeof(double));
    integer *iwork = calloc(n,sizeof(integer));
    dtrcon_("1","L","N",&n,A,&n,&rcond,work,iwork,&info);
    printf("%lf\n",1/rcond);
    dtrcon_("I","L","N",&n,A,&n,&rcond,work,iwork,&info);
    printf("%lf\n",1/rcond);

    printf("System 2\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A_0[i+j*n] = A_0[i+j*n]*1e-06;
        }
    }
    vector_copy(n*n,A_0,A);
    vector_copy(n,b,x);
    print_matrix(n,A);
    print_vector(n,x);

    printf("solving system:\n");

    printf("solving with dtrtrs:\n");
    dtrtrs_("L","N","N", &n, &nhrs, A, &n, x, &n, &info);
    print_vector(n,x);
    printf("solving with dgesv:\n");
    vector_copy(n,b,x);
    dgesv_(&n, &nhrs, A, &n, ipvt, x, &n, &info);
    print_vector(n,x);

    vector_copy(n*n,A_0,A);

    printf("Condition numbers\n");
    work = calloc(5*n, sizeof(double));
    dgesvd_("N","N", &n, &n, A, &n, s, NULL, &n, NULL, &n, work, &lwork, &info);
    printf("%lf\n", s[0]/s[n-1]); //compute 2-norm condition number
    print_vector(n,s);

    vector_copy(n*n,A_0,A);
    free(work);
    work = calloc(3*n,sizeof(double));
    dtrcon_("1","L","N",&n,A,&n,&rcond,work,iwork,&info);
    printf("%lf\n",1/rcond);
    dtrcon_("I","L","N",&n,A,&n,&rcond,work,iwork,&info);
    printf("%lf\n",1/rcond);


    return 0;
}

void print_vector(integer n, doublereal *x)
{
    integer i;
    for (i = 0; i < n; i++) {
        printf("%f ", x[i]);
    }
    printf("\n\n");
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
