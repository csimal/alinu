#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 * prints a vector
 * IN:
 * n : the size of the vector
 * x : a pointer to a bloc of size n of double containing the elements of the vector
 * */
void print_vector(int n, double *x)
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%le ", x[i]);
    }
    printf("\n\n");
}

/*
 * prints a matrix A
 * IN:
 * m : the number of rows of the matrix
 * n : the number of columns of the matrix
 * a : a pointer to a bloc of size m*n of double such that
 *      A[i][j] = a[j*m +i], (0<=j<n, 0<=i<m)
 * */
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

/*
 * read a vector from a file
 * IN:
 * fp : a pointer to an opened file where the vector to be read is on the current line
 * n : the number of elements of the vector
 * v : a pointer to a bloc of size n of doubles. On exit, contains the vector read from the file
 * */
void read_vector(FILE* fp, int n, double *v) {
    char line[1024], *p, *e;
    int i;

    fgets(line, sizeof(line), fp);
    for (i = 0, p = line; i < n; i++, p = e) {
        v[i] = strtod(p, &e);
    }
}

/*
 * read a matrix A from a file and format it in column major order
 * IN:
 * fp : a pointer to an opened file such that the m next lines of the file contain the m lines of A
 * m : the number of rows of A
 * n : the number of columns of A
 * a : a pointer to a bloc of size m*n of doubles. On exit, contains the matrix A in column major order
 *      i.e. A[i][j] = A[i+j*m], (0<=i<m, 0<=j<n)
 * */
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

/* 
 * read a symmetric (or lower triangular) matrix A in column form
 * IN:
 * fp : a pointer to an opened file such that the successive elements of A
 *      are stored column first on each line.
 * n : the number of rows and columns of A
 * a_vec : a pointer to a bloc of size (n*(n+1))/2 of doubles.
 *      On exit, contains A in packed storage, i.e.
 *      A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0 <= j <= i < n)
 * */
void read_matrix_sc(FILE* fp, int n, double *a_vec) {
    int i,j;
    double val;

    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            fscanf(fp,"%lf\n", &val);
            if (j <= i) {
                a_vec[j*n + i - (j*(j+1))/2] = val;
            }
        }
    }
}

/*
 * write a vector in column form to a file
 * IN:
 * fp : a pointer to an opened file where the vector should be written
 * n : the number if elements of the vector
 * v : a pointer to a bloc of size n of doubles containing the elements of the vector
 * */
void write_vector_c(FILE *fp, int n, double *v) {
    int i;
    fprintf(fp, "%d\n",n);
    for (i = 0; i < n; i++) {
       fprintf(fp, "%le\n",v[i]); 
    }
}

/*
 * write a vector in row form to a file
 * IN:
 * fp : a pointer to an opened file where the vector should be written
 * n : the number if elements of the vector
 * v : a pointer to a bloc of size n of doubles containing the elements of the vector
 * */
void write_vector(FILE *fp, int n, double *v) {
    int i;
    for (i = 0; i < n; i++) {
        fprintf(fp, "%le ", v[i]);
    }
    fprintf(fp, "\n");
}

/*
 * copy a vector to another
 * IN:
 * n : the size of the vectors
 * u : a pointer to a bloc of size n of doubles containing the vector to be copied
 * v : a pointer to a bloc of size n of doubles. On exit, the elements of u are copied to v
 *      i.e. v[i] = u[i], (0<=i<n)
 * */
void copy_vector(int n, double *u, double *v) {
    int i;
    for (i = 0; i < n; i++) {
        v[i] = u[i];
    }
}

/*
 * computes A*x + y, where A is a symmetric matrix in packed storage
 * IN:
 * n : the number of rows and columns of A
 * a_vec : a pointer to a bloc of size (n*(n+1))/2 of doubles containing A in packed storage.
 *      i.e. A[i][j] = a_vec[j*n + i - (j*(j+1))/2], (0<=j<=i<n)
 * x : a pointer to a bloc of size n of doubles containing the elements of x
 * y : a pointer to a bloc of size n of doubles containing the elements of y
 *      On exit, contains the result of A*x + y
 * */
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

/*
 * compute a*x + y, where a is a scalar
 * IN:
 * n : the number of elements of x and y
 * a : a real number
 * x : a pointer to a bloc of size n of doubles containing the elements of x
 * y : a pointer to a bloc of size n of doubles containing the elements of y
 *      On exit, contains the result of a*x + y
 * */
void saxpy(int n, double a, double *x, double *y) {
    int i;
    for (i = 0; i < n; i++) {
        y[i] += a*x[i];
    }
}

/*
 * set all the elements of a vector to 0
 * IN:
 * n : the number of elements of the vector
 * v : a pointer to a bloc of size n of doubles. On exit, v[i] = 0.0 (0<=i<n)
 * */
void vector_nullify(int n, double *v) {
    int i;
    for (i = 0; i < n; i++) {
        v[i] = 0.0;
    }
}

/*
 * compute the dot product of u and v using Kahan sum algorithm
 * IN:
 *  n : the number of elements of u and v
 *  u : a pointer to a bloc of size n of doubles containing the elements of u
 *  v : a pointer to a bloc of size n of doubles containing the elements of v
 *  OUT:
 *  returns the dot product u[0]*v[0] + ... + u[n-1]*v[n-1]
 */
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

/*
 * compute the euclidian norm of a vector
 * IN:
 * n : the number of elements of the vector
 * v : a pointer to a bloc of size n of doubles containing the elements of the vector
 * OUT:
 * returns the euclidian norm of v
 * */
double vector_norm(int n, double *v) {
    return sqrt(vector_dot_product(n,v,v));
}

/*
 * Converts a symmetric matrix from packed storage to full storage
 * IN:
 * n : the order the matrix
 * a_vec : the matrix in packed storage, i.e. a pointer to a bloc of size (n*(n+1))/2 of doubles
 *      such that A[i][j] = a_vec[j*n + i -(j*(j+1))/2], (0 <= j <= i < n)
 * OUT:
 * returns a pointer to a bloc of size n*n of doubles such that
 * A[i][j] = a_full[j*n+i], 0 <= i < n, 0 <= j < n
 * */
double* packed_to_full(int n, double *a_vec) {
    double *a_full = calloc(n*n, sizeof(double));
    int i,j;

    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            a_full[j*n+i] = a_vec[j*n + i - (j*(j+1))/2];
            a_full[i*n+j] = a_vec[j*n + i - (j*(j+1))/2];
        }
    }
    return a_full;
}

/*
 * Converts a symmetric matrix from full storage to packed storage
 * IN:
 * n : the order of the matrix
 * a_full : a pointer to a bloc of size n*n of doubles such that
 *       A[i][j] = a_full[j*n+i], 0 <= i < n, 0 <= j < n
 * OUT:
 * returns a pointer to a bloc of size (n*(n+1))/2 of doubles
 * such that A[i][j] = a_vec[j*n + i -(j*(j+1))/2], (0 <= j <= i < n)
 * */
double* full_to_packed(int n, double *a_full) {
    double *a_vec = calloc((n*(n+1))/2,sizeof(double));
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            a_vec[j*n + i - (j*(j+1))/2] = a_full[j*n+i];
        }
    }
    return a_vec;
}
