#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "matrix.h"

#define TOL 0.001

void solve_cg(int n, double *a, double *b, double *x, double tol);
void solve_pcg(int n, double *a, double *b, double *m_vec);

int main(int argc, char const *argv[]) {
    FILE* fp;
    double *a_vec, *b, *x;
    int i,n;
    char filename[16];
    int lcm = LAPACK_COL_MAJOR;

    fp = fopen("donneeRHS.dat","r");
    fscanf(fp,"%d\n",&n);

    b = calloc(n,sizeof(double));
    read_matrix(fp,n,1,b);
    print_vector(n,b);
    a_vec = calloc(n*(n+1)/2,sizeof(double));
    x = calloc(n,sizeof(double));

    fclose(fp);

    for (i = 1; i < 4; i++) {
        sprintf(filename,"donneeSyst%d.dat",i);
        fp = fopen(filename,"r");
        fscanf(fp, "%d\n", &n);
        read_matrix_sc(fp, n, a_vec);
        solve_cg(n, a_vec, b, x, TOL);
        print_vector(n,x);
        double *y = calloc(n,sizeof(double));
        gaxpy_s(n, a_vec, x, y);
        saxpy(n, -1.0, b, y);
        printf("norm(A*x-b) = %lf\n", vector_norm(n,y));
        vector_nullify(n,x);
    }


    return 0;
}


void solve_cg(int n, double *a_vec, double *b, double *x, double tol) {
    double *r = calloc(n,sizeof(double));
    double *p = calloc(n,sizeof(double));
    int i,j,k;

    saxpy(n, -1.0, b, r); // r <- -b
    gaxpy_s(n, a_vec, x, r); // r <- A*x -b
    saxpy(n, -1.0, r, p); // p <- -r

    double alpha, *ap, pap, beta;
    ap = calloc(n,sizeof(double));

    k = 0;
    while (vector_norm(n,r) >= tol && k++<n) {
        printf("k = %d\n",k-1);
        gaxpy_s(n, a_vec, p, ap); // compute A*p_k only once
        pap = vector_dot_product(n, p, ap); // same with p_k'*A*p_k
        printf("pap = %lf\n",pap);
        alpha = -vector_dot_product(n,r,r)/pap;
        printf("alpha = %lf\n",alpha);
        saxpy(n, alpha, p, x); // x_k+1 = x_k+1 + alpha*p_k
        saxpy(n, -alpha, ap, r); // r_k+1 = r_k - alpha*A*p_k
        beta = vector_dot_product(n, r, ap)/pap;

        for (i = 0; i < n; i++) {
            p[i] = beta*p[i] - r[i];
        }

        vector_nullify(n,ap); // set ap to zero for next iteration
    }
    free(r);
    free(p);
    free(ap);
}

void solve_pcg(int n, double *a_vec, double *b, double *m_vec) {

}
