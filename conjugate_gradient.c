#include "conjugate_gradient.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "matrix.h"

#define LCM LAPACK_COL_MAJOR

/*
 * Solve a symmetric definite positive linear system A*x = b using the conjugate gradient
 * method without preconditionner. Writes the iterations to a given file
 * IN:
 * fp : a pointer to an opened file where each iterate x_k is to be written
 * n : the order of the system
 * a_vec : the matrix A in packed storage. i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * b : a pointer to a bloc of size n of doubles containing the right hand side of the equation
 * x : a pointer to a bloc of size n of doubles containing an arbitrary initial guess of the solution
 * tol : the tolerance at which to stop the algorithm. 0 < tol
 * OUT:
 * returns the number of iterations at which the algorithm stopped
 * */
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

/*
 * Solve a symmetric definite positive linear system A*x = b using the conjugate gradient
 * method with the jacobi preconditionner. Writes the iterations to a given file
 * The jacobi preconditionner is the diagonal matrix whose elements are those of the diagonal of A.
 * IN:
 * fp : a pointer to an opened file where each iterate x_k is to be written
 * n : the order of the system
 * a_vec : the matrix A in packed storage. i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * b : a pointer to a bloc of size n of doubles containing the right hand side of the equation
 * x : a pointer to a bloc of size n of doubles containing an arbitrary initial guess of the solution
 * tol : the tolerance at which to stop the algorithm. 0 < tol
 * OUT:
 * returns the number of iterations at which the algorithm stopped
 * */
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

/*
 * Solve a symmetric definite positive linear system A*x = b using the conjugate gradient
 * method with the SSOR preconditionner. Writes the iterations to a given file
 * The SSOR preconditionner is defined as M(w) = (1/(2-w))*(D/w + L)*(D/w)^-1 *(D/w + L)'
 * where A = L + D + L', D is diagonal, L is a lower triangular matrix.
 * IN:
 * fp : a pointer to an opened file where each iterate x_k is to be written
 * n : the order of the system
 * a_vec : the matrix A in packed storage. i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * b : a pointer to a bloc of size n of doubles containing the right hand side of the equation
 * x : a pointer to a bloc of size n of doubles containing an arbitrary initial guess of the solution
 * w : the parameter used in the SSOR preconditionner. 0 < w < 2
 * tol : the tolerance at which to stop the algorithm. 0 < tol
 * OUT:
 * returns the number of iterations at which the algorithm stopped
 * */
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
    for (i = 0; i < n; i++) {
        y[i] *= a_vec[i*n + i - (i*(i+1))/2]/((2.0-w)*w); // y = (D/((2-w)*w))*s
    }
    LAPACKE_dtptrs(LCM, 'L', 'T', 'N', n, 1, m_vec, y, n); // M'*s = q

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
            y[i] *= a_vec[i*n + i - (i*(i+1))/2]/((2.0-w)*w); // y = (D/((2-w)*w))*s
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

/*
 * Solve a symmetric definite positive linear system A*x = b using the conjugate gradient
 * method with the incomplete cholesky preconditionner. Writes the iterations to a given file
 * The incomplete cholesky preconditionner is defined as M = K*K'
 * where K is a lower triangular matrix such that K*K' ~= A. i.e. K is an incomplete cholesky factorisation of A
 * IN:
 * fp : a pointer to an opened file where each iterate x_k is to be written
 * n : the order of the system
 * a_vec : the matrix A in packed storage. i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * b : a pointer to a bloc of size n of doubles containing the right hand side of the equation
 * x : a pointer to a bloc of size n of doubles containing an arbitrary initial guess of the solution
 * k_vec : an incomplete cholesky factorization of A in packed storage. i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that K[i][j] = k_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * tol : the tolerance at which to stop the algorithm. 0 < tol
 * OUT:
 * returns the number of iterations at which the algorithm stopped
 * */
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

/*
 * Solve a symmetric definite positive linear system A*x = b using the conjugate gradient
 * method with the spectral preconditionner of order m. Writes the iterations to a given file
 * The spectral conditionner of order m is defined as M = I + (lambda[1]-1)*v[1]*v[1]' + ... + (lambda[m]-1)*v[m]*v[m]'
 * where lambda[i] is the i-th smallest eigenvalue of A, and v[i] is the associated eigenvector.
 * IN:
 * fp : a pointer to an opened file where each iterate x_k is to be written
 * n : the order of the system
 * a_vec : the matrix A in packed storage. i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * b : a pointer to a bloc of size n of doubles containing the right hand side of the equation
 * x : a pointer to a bloc of size n of doubles containing an arbitrary initial guess of the solution
 * m : the order of the preconditionner
 * tol : the tolerance at which to stop the algorithm. 0 < tol
 * OUT:
 * returns the number of iterations at which the algorithm stopped
 * */
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

/*
 * Computes the m smallest eigenvalues of a symmetric matrix A, and optionnally the associated eigenvectors
 * IN:
 * n : the order of A
 * a_vec : the matrix A in packed storage, i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * m : the desired number of eigenvalues
 * jobz : Parameter for computing eigenvectors as well. 'N' for eigenvalues only, 'V' for eigenvalues and eigenvectors
 * v : a pointer to a bloc of size n*m of doubles. If jobz = 'V', contains the computed eigenvectors.
 *      if jobz = 'N', it will not be referenced.
 * OUT:
 * returns a pointer to a bloc of size n of double such that it's i-th element contains the i-th smallest eigenvalue
 * of A. (0 <= i < m)
 * */
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

void write_eig_jacobi(FILE *fp, int n, double *a_vec) {
    double *t_vec = calloc((n*(n+1))/2,sizeof(double));
    int i,j;

    copy_vector((n*(n+1))/2, a_vec, t_vec);

    for (i = 0; i < n; i++) {
        double d = a_vec[i*n + i - (i*(i+1))/2];
        for (j = 0; j <= i; j++) {
            t_vec[j*n + i - (j*(j+1))/2] /= d; // divide each line by it's diagonal element
        }
    }

    double *lambda = eig(n, t_vec, n, 'N', NULL);
    write_vector_c(fp, n, lambda);

    free(t_vec);
    free(lambda);
}

void write_eig_ssor(FILE *fp, int n, double *a_vec, double w) {
    double *t_vec = calloc((n*(n+1))/2,sizeof(double));
    int i,j;

    copy_vector((n*(n+1))/2, a_vec, t_vec);

    for (i = 0; i < n; i++) {
        t_vec[i*n + i - (i*(i+1))/2] /= w;
    }

    double *a_full = packed_to_full(n, a_vec);

    LAPACKE_dtptrs(LCM, 'L', 'T', 'N', n, n, t_vec, a_full, n);
    for (i = 0; i < n; i++) {
        double d = a_full[i*n+i]/((2.0-w)*w);
        for (j = 0; j < n; j++) {
            a_full[j*n+i] *= d;
        }
    }
    LAPACKE_dtptrs(LCM, 'L', 'N', 'N', n, n, t_vec, a_full, n);
    double *m_vec = full_to_packed(n, a_full);

    double *lambda = eig(n, m_vec, n, 'N', NULL);
    write_vector_c(fp, n, lambda);

    free(lambda);
    free(m_vec);
    free(a_full);
    free(t_vec);
}

void write_eig_cholesky(FILE *fp, int n, double *a_vec, double *k_vec) {
    double *a_full = packed_to_full(n, a_vec);

    LAPACKE_dtptrs(LCM, 'L', 'T', 'N', n, n, k_vec, a_full, n);
    LAPACKE_dtptrs(LCM, 'L', 'N', 'N', n, n, k_vec, a_full, n);

    double *m_vec = full_to_packed(n, a_full);

    double *lambda = eig(n, m_vec, n, 'N', NULL);
    write_vector_c(fp, n, lambda);

    free(lambda);
    free(m_vec);
    free(a_full);
}

void write_eig_spectral(FILE *fp, int n, double *a_vec, int m) {
    double *m_vec = calloc((n*(n+1))/2,sizeof(double));
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
    free(v);

    double *a_full = packed_to_full(n, a_vec);

    LAPACKE_dppsv(LCM, 'L', n, n, m_vec, a_full, n);

    free(m_vec);
    m_vec = full_to_packed(n, a_full);

    lambda = eig(n, m_vec, n, 'N', NULL);
    write_vector_c(fp, n, lambda);

    free(lambda);
    free(m_vec);
    free(a_full);
}
