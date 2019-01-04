#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "matrix.h"
#include "conjugate_gradient.h"

#define TOL 0.001
#define LCM LAPACK_COL_MAJOR

int main(int argc, char const *argv[]) {
    FILE *fp, *fp2;
    double *a_vec, *a_tmp, *b, *x, *lambda;
    int i,n,k;
    char filename[100];

    fp = fopen("donneeRHS.dat","r");
    fscanf(fp,"%d\n",&n);
    b = calloc(n,sizeof(double));
    read_matrix(fp,n,1,b);
    fclose(fp);

    a_vec = calloc((n*(n+1))/2,sizeof(double));
    a_tmp = calloc((n*(n+1))/2,sizeof(double));
    x = calloc(n,sizeof(double));

    for (i = 1; i < 4; i++) {
        // read A
        printf("Processing matrix A_%d\n",i);
        sprintf(filename,"donneeSyst%d.dat",i);
        fp = fopen(filename,"r");
        fscanf(fp, "%d\n", &n);
        read_matrix_sc(fp, n, a_vec);
        fclose(fp);

        // compute and save eigenvalues of A
        printf("Computing eigenvalues\n");
        copy_vector((n*(n+1))/2, a_vec, a_tmp);
        lambda = eig(n, a_tmp, n, 'N', NULL); // it's ok to pass NULL as it won't be referenced
        sprintf(filename, "results/eigenvalues_A_%d.dat", i);
        fp = fopen(filename,"w");
        write_vector_c(fp, n, lambda);
        free(lambda);
        fclose(fp);
        // solving with direct method (cholesky) for comparison
        printf("Computing direct solution\n");
        copy_vector(n,b,x);
        char equed;
        double rcond, ferr, berr;
        LAPACKE_dppsvx(LCM, 'N', 'L', n, 1, a_vec, a_tmp, &equed, NULL, b, n, x, n, &rcond, &ferr, &berr);
        sprintf(filename,"results/direct_solution_%d.dat", i);
        fp = fopen(filename,"w");
        write_vector_c(fp, n, x);
        fclose(fp);
        printf("Estimate of K(A_%d): %lf\n", i, rcond);
        printf("Error bound on solution: %lf\n", ferr);
        
        // solve with non-preconditionned cg
        printf("Solving with non-preconditionned cg\n");
        copy_vector(n,b,x);
        sprintf(filename, "results/termes_cg_%d.dat",i);
        fp = fopen(filename,"w");
        k = solve_cg(fp, n, a_vec, b, x, TOL);
        fclose(fp);
        printf("converged after %d iterations\n", k);

        // solve with cg with jacobi precondtionner
        printf("Solving with Jacobi preconditionner\n");
        copy_vector(n,b,x);
        sprintf(filename, "results/termes_cg_jacobi_%d.dat", i);
        fp = fopen(filename,"w");
        k = solve_cg_jacobi(fp, n, a_vec, b, x, TOL);
        fclose(fp);
        printf("converged after %d iterations\n", k);

        sprintf(filename, "results/eigenvalues_jacobi_%d.dat", i);
        fp = fopen(filename,"w");
        write_eig_jacobi(fp, n, a_vec);
        fclose(fp);
        
        // solve with cg with ssor precondtionner
        printf("Solving with SSOR preconditionner\n");
        double w;
        sprintf(filename, "results/iter_ssor_%d.dat", i);
        fp2 = fopen(filename, "w");
        for (w = 0.1; w < 2.0; w += 0.1) {
            copy_vector(n,b,x);
            sprintf(filename, "results/termes_cg_ssor_w%.2lf_%d.dat", w, i);
            fp = fopen(filename,"w");
            k = solve_cg_ssor(fp, n, a_vec, b, x, w, TOL);
            fclose(fp);
            printf("converged after %d iterations\n", k);
            fprintf(fp2, "%.2lf %d\n", w, k);

            sprintf(filename, "results/eigenvalues_ssor_w%.2lf_%d.dat", w, i);
            fp = fopen(filename,"w");
            write_eig_ssor(fp, n, a_vec, w);
            fclose(fp);
        }
        fclose(fp2);

        // solve with cg with incomplete cholesky precondtionner
        printf("Solving with incomplete Cholesky conditionner\n");
        sprintf(filename, "donneeChol%d.dat", i);
        fp = fopen(filename, "r");
        read_matrix_sc(fp, n, a_tmp); // a_tmp stands in for k_vec
        fclose(fp);
        copy_vector(n,b,x);
        sprintf(filename, "results/termes_cg_ichol_%d.dat", i);
        fp = fopen(filename,"w");
        k = solve_cg_cholesky(fp, n, a_vec, b, x, a_tmp, TOL);
        fclose(fp);
        printf("converged after %d iterations\n", k);
        
        sprintf(filename, "results/eigenvalues_cholesky_%d.dat", i);
        fp = fopen(filename,"w");
        write_eig_cholesky(fp, n, a_vec, a_tmp);
        fclose(fp);
        
        // solve with cg with spectral precondtionner
        printf("Solving with spectral conditionner\n");
        sprintf(filename, "results/iter_spectral_%d.dat", i);
        fp2 = fopen(filename, "w");
        int m;
        for (m = 1; m <= 100; m++) {
            printf("m = %d\n", m);
            copy_vector(n,b,x);
            sprintf(filename, "results/termes_cg_spectral_%d_%d.dat", m, i);
            fp = fopen(filename,"w");
            k = solve_cg_spectral(fp, n, a_vec, b, x, m, TOL);
            fclose(fp);
            printf("converged after %d iterations\n", k);
            fprintf(fp2, "%d %d\n", m, k);
            sprintf(filename, "results/eigenvalues_spectral_%d_%d.dat", m, i);
            fp = fopen(filename,"w");
            write_eig_spectral(fp, n, a_vec, m);
            fclose(fp);
        }
        fclose(fp2);
    }

    free(a_vec);
    free(a_tmp);
    free(b);
    free(x);

    return 0;
}
