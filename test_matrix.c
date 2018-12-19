#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main(int argc, char const *argv[]) {
    FILE *fp;
    char *mode = "r";
    fp = fopen("donnees_Lxb.txt",mode);
    integer n;
    doublereal *a;
    fscanf(fp,"%ld", &n);
    a = calloc(n*n,sizeof(doublereal));
    read_matrix(fp,n,n,a);
    print_matrix(n,a);
    return 0;
}
