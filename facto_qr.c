#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "cblas.h"
#include "matrix.h"

int main(int argc, const char* argv[])
{
    FILE *fp;
    int m,n,i,j,k,l,nhrs=1;
    doublereal *A, *x, *b, *tau, *work, *Q;
    int info, lwork;
    char filename[12];

     
    for (l = 0; l < 3; l++) {
        sprintf(filename, "donnee%d.dat",l+1);
        fp = fopen(filename,"r");
        fscanf(fp,"%d\n",&m);
        fscanf(fp,"%d\n",&n);
        printf("m=%d, n=%d\n",m,n);

        lwork = n;

        A = calloc(m*n,sizeof(doublereal));
        tau = calloc(n,sizeof(doublereal));
        work = calloc(lwork,sizeof(doublereal));
        b = calloc(m,sizeof(doublereal));
        read_matrix(fp,m,n,A);
        read_matrix(fp,m,1,b);

        printf("A%d\n",l+1);
        print_matrix(m,n,A);
        printf("b%d\n",l+1);
        print_vector(m,b);
        
        dgeqrf_(&m, &n, A, &m, tau, work, &lwork, &info);
        printf("optimal lwork: %lf\n",work[0]);
        Q = (doublereal *) extraireQ(m, n, A, tau);
        printf("Q%d:\n",l+1);
        print_matrix(m,m,Q);

        x = calloc(n,sizeof(doublereal));
        // compute Q'*b[1:n]
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                x[i] += Q[j+i*m]*b[j];
            }
        }
        printf("c%d:\n",l+1);
        print_vector(n,x);
        dtrtrs_("U","N","N",&n,&nhrs,A,&m,x,&n,&info);
        printf("x%d:\n",l+1);
        print_vector(n,x);


        fclose(fp);
    }
   
    return 0;
}

