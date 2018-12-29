#include "conjugate_gradient.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "matrix.h"

#define LCM LAPACK_COL_MAJOR

int solve_cg(FILE *fp, int n, double *a_vec, double *b, double *x, double tol) {
    double *r = calloc(n,sizeof(double));
    double *p = calloc(n,sizeof(double));
    int i,k;

    write_vector(fp, n, x);

    gaxpy_s(n, a_vec, x, r); // r <- A*x
    saxpy(n, -1.0, b, r); // r <- A*x-b
    for (i = 0; i < n; i++) {
        r[i] = -r[i]; // r <- b-A*x
    }
    copy_vector(n,r,p); // p <- r

    double alpha, *ap, pap, beta;
    ap = calloc(n,sizeof(double));

    k = 0;
    while (vector_norm(n,r) >= tol) {
        gaxpy_s(n, a_vec, p, ap); // compute A*p_k only once
        pap = vector_dot_product(n, p, ap); // same with p_k'*A*p_k
        alpha = vector_dot_product(n,r,r)/pap;
        saxpy(n, alpha, p, x); // x_k+1 = x_k + alpha*p_k
        saxpy(n, -alpha, ap, r); // r_k+1 = r_k - alpha*A*p_k
        beta = vector_dot_product(n, r, ap)/pap;

        for (i = 0; i < n; i++) {
            p[i] = beta*p[i] + r[i];
        }

        vector_nullify(n,ap); // set ap to zero for next iteration
        ++k;
        write_vector(fp, n, x);
    }
    free(r);
    free(p);
    free(ap);

    return k;
}

int solve_cg_jacobi(FILE *fp, int n, double *a_vec, double *b, double *x, double tol) {
    double *r = calloc(n,sizeof(double));
    double *y = calloc(n,sizeof(double));
    double *p = calloc(n,sizeof(double));
    int i,k;

    write_vector(fp, n, x);

    gaxpy_s(n, a_vec, x, r); // r <- A*x
    saxpy(n, -1.0, b, r); // r <- A*x-b
    for (i = 0; i < n; i++) {
        r[i] = -r[i]; // r <- b-A*x
        y[i] = r[i]/a_vec[i*n+i-(i*(i+1))/2]; // M*y = r
        p[i] = y[i];
    }

    double alpha, *ap, pap, beta, yr;
    ap = calloc(n,sizeof(double));
    k = 0;
    while (vector_norm(n,r) >= tol) {
        gaxpy_s(n, a_vec, p, ap); // ap <- A*p_k
        pap = vector_dot_product(n, p, ap); // pap <- p_k'*A*p_k
        yr = vector_dot_product(n,y,r); // yr <- y_k'*r_k
        alpha = yr/pap;
        saxpy(n, alpha, p, x); // x_k+1 = x_k + alpha*p_k
        saxpy(n, -alpha, ap, r); // r_k+1 = r_k - alpha*A*p_k

        for (i = 0; i < n; i++) {
            y[i] = r[i]/a_vec[i*n+i-(i*(i+1))/2]; // M*y = r
        }

        beta = vector_dot_product(n, y, r)/yr;
        
        for (i = 0; i < n; i++) {
            p[i] = y[i] + beta*p[i];
        }
        
        vector_nullify(n,ap);
        write_vector(fp,n,x);
        ++k;
    }
    free(r);
    free(y);
    free(p);
    free(ap);

    return k;
}

int solve_cg_ssor(FILE *fp, int n, double *a_vec, double *b, double *x, double w, double tol) {
    double *r = calloc(n,sizeof(double));
    double *y = calloc(n,sizeof(double));
    double *p = calloc(n,sizeof(double));
    double *m_vec = calloc((n*(n+1))/2,sizeof(double));
    int i,k;

    write_vector(fp, n, x);

    copy_vector((n*(n+1))/2, a_vec, m_vec);
    for (i = 0; i < n; i++) {
        m_vec[i*n + i - (i*(i+1))/2] /= w; // M = (D/w + L)
    }

    gaxpy_s(n, a_vec, x, r); // r <- A*x
    saxpy(n, -1.0, b, r); // r <- A*x-b
    for (i = 0; i < n; i++) {
        r[i] = -r[i]; // r_0 =  b-A*x_0
    }

    copy_vector(n,r,y);
    LAPACKE_dtptrs(LCM, 'L', 'N', 'N', n, 1, m_vec, y, n); // M*q = r
    LAPACKE_dtptrs(LCM, 'L', 'T', 'N', n, 1, m_vec, y, n); // M'*s = q
    for (i = 0; i < n; i++) {
        y[i] *= w*(2.0-w)*a_vec[i*n + i - (i*(i+1))/2]; // y = w*(2-w)*D*s
    }

    copy_vector(n,y,p); // p_0 = y_0

    double alpha, *ap, pap, beta, yr;
    ap = calloc(n,sizeof(double));
    k = 0;
    while (vector_norm(n,r) >= tol) {
        gaxpy_s(n, a_vec, p, ap); // ap <- A*p_k
        pap = vector_dot_product(n, p, ap); // pap <- p_k'*A*p_k
        yr = vector_dot_product(n,y,r); // yr <- y_k'*r_k
        alpha = yr/pap;
        saxpy(n, alpha, p, x); // x_k+1 = x_k + alpha*p_k
        saxpy(n, -alpha, ap, r); // r_k+1 = r_k - alpha*A*p_k

        copy_vector(n,r,y);
        LAPACKE_dtptrs(LCM, 'L', 'T', 'N', n, 1, m_vec, y, n); // M'*q = r
        for (i = 0; i < n; i++) {
            y[i] *= w*(2-w)*a_vec[i*n + i - (i*(i+1))/2]; // (D/w)^-1*s = q
        }
        LAPACKE_dtptrs(LCM, 'L', 'N', 'N', n, 1, m_vec, y, n); // M*y = s
        beta = vector_dot_product(n, y, r)/yr;

        for (i = 0; i < n; i++) {
            p[i] = y[i] + beta*p[i];
        }

        vector_nullify(n,ap);
        write_vector(fp,n,x);
        ++k;
    }
    free(r);
    free(y);
    free(p);
    free(ap);
    free(m_vec);

    return k;
}

int solve_cg_cholesky(FILE *fp, int n, double *a_vec, double *b, double *x, double *k_vec, double tol) {
    double *r = calloc(n,sizeof(double));
    double *y = calloc(n,sizeof(double));
    double *p = calloc(n,sizeof(double));
    int i,k;

    write_vector(fp, n, x);

    gaxpy_s(n, a_vec, x, r); // r <- A*x
    saxpy(n, -1.0, b, r); // r <- A*x-b
    for (i = 0; i < n; i++) {
        r[i] = -r[i]; // r_0 =  b-A*x_0
    }

    copy_vector(n,r,y);
    LAPACKE_dtptrs(LCM, 'L', 'T', 'N', n, 1, k_vec, y, n); // K'*q = r
    LAPACKE_dtptrs(LCM, 'L', 'N', 'N', n, 1, k_vec, y, n); // K*y = q

    copy_vector(n,y,p); // p_0 = y_0

    double alpha, *ap, pap, beta, yr;
    ap = calloc(n,sizeof(double));
    k = 0;
    while (vector_norm(n,r) >= tol) {
        gaxpy_s(n, a_vec, p, ap); // ap <- A*p_k
        pap = vector_dot_product(n, p, ap); // pap <- p_k'*A*p_k
        yr = vector_dot_product(n,y,r); // yr <- y_k'*r_k
        alpha = yr/pap;
        saxpy(n, alpha, p, x); // x_k+1 = x_k + alpha*p_k
        saxpy(n, -alpha, ap, r); // r_k+1 = r_k - alpha*A*p_k

        copy_vector(n,r,y);
        LAPACKE_dtptrs(LCM, 'L', 'T', 'N', n, 1, k_vec, y, n); // K'*q = r
        LAPACKE_dtptrs(LCM, 'L', 'N', 'N', n, 1, k_vec, y, n); // K*s = q

        beta = vector_dot_product(n, y, r)/yr;

        for (i = 0; i < n; i++) {
            p[i] = y[i] + beta*p[i];
        }

        vector_nullify(n,ap);
        write_vector(fp,n,x);
        ++k;
    }
    free(r);
    free(y);
    free(p);
    free(ap);

    return k;
}

int solve_cg_spectral(FILE *fp, int n, double *a_vec, double *b, double *x, int m, double tol) {
    double *r = calloc(n,sizeof(double));
    double *y = calloc(n,sizeof(double));
    double *p = calloc(n,sizeof(double));
    double *m_vec = calloc((n*(n+1))/2,sizeof(double));
    double *t_vec = calloc((n*(n+1))/2,sizeof(double));
    double *v = calloc(n*m, sizeof(double)); // eigenvectors
    double *lambda;
    int i,j,k;
    
    lambda = eig(n, a_vec, m, 'V', v); // get m smallest eigenvalues and eigenvectors

    for (j = 0; j < n; j++) {
        m_vec[j*n + j - (j*(j+1))/2] = 1.0;
        for (i = j; i < n; i++) {
            for (k = 0; k < m; k++) {
                m_vec[j*n + i - (j*(j+1))/2] += (lambda[k]-1)*v[k*n+i]*v[k*n+j];
            }
        }
    }
    free(lambda);
    free(v); // won't be needing those anymore

    write_vector(fp, n, x);

    gaxpy_s(n, a_vec, x, r); // r <- A*x
    saxpy(n, -1.0, b, r); // r <- A*x-b
    for (i = 0; i < n; i++) {
        r[i] = -r[i]; // r_0 =  b-A*x_0
    }

    copy_vector(n,r,y);
    copy_vector((n*(n+1))/2, m_vec, t_vec); // must copy because dppsv modifies it's arguments
    LAPACKE_dppsv(LCM, 'L', n, 1, t_vec, y, n); // M*y = r

    copy_vector(n,y,p); // p_0 = y_0

    double alpha, *ap, pap, beta, yr;
    ap = calloc(n,sizeof(double));
    k = 0;
    while (vector_norm(n,r) >= tol) {
        gaxpy_s(n, a_vec, p, ap); // ap <- A*p_k
        pap = vector_dot_product(n, p, ap); // pap <- p_k'*A*p_k
        yr = vector_dot_product(n,y,r); // yr <- y_k'*r_k
        alpha = yr/pap;
        saxpy(n, alpha, p, x); // x_k+1 = x_k + alpha*p_k
        saxpy(n, -alpha, ap, r); // r_k+1 = r_k - alpha*A*p_k

        copy_vector(n,r,y);
        copy_vector((n*(n+1))/2, m_vec, t_vec); // must copy because dppsv modifies it's arguments
        LAPACKE_dppsv(LCM, 'L', n, 1, t_vec, y, n); // M*y = r

        beta = vector_dot_product(n, y, r)/yr;

        for (i = 0; i < n; i++) {
            p[i] = y[i] + beta*p[i];
        }

        vector_nullify(n,ap);
        write_vector(fp,n,x);
        ++k;
    }
    free(r);
    free(y);
    free(p);
    free(ap);
    free(m_vec);
    free(t_vec);

    return k;
}

double* eig(int n, double *a_vec, int m, char jobz, double *v) {
    double *lambda = calloc(n,sizeof(double));
    double *t_vec = calloc((n*(n+1))/2, sizeof(double));
    int k;
    int *ifail = calloc(n,sizeof(int));
    double abstol = LAPACKE_dlamch_work('S'); // optimal tolerance parameter

    copy_vector((n*(n+1))/2, a_vec, t_vec); // must copy because dspevx modifies it

    LAPACKE_dspevx(LCM, jobz, 'I', 'L', n, t_vec, 0.0, 0.0, 1, m, abstol, &k, lambda, v, n, ifail);
    
    free(t_vec);
    free(ifail);
    return lambda;
}
