#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cblas.h>
#include "f2c.h"

/* Fonction extraireQ(int m, int n, double *A, double *tau)
	Auteur : C. Beauthier (Dec. 2008) - Code modifie par J. Blanchard (Nov. 2017)
	But : Fonction qui permet d'extraire de la matrice A (mxn)(donnee par la fonction DGEQRF)
		la matrice Q produit de reflecteurs
	 Q = H1 H2... Hn ou Hi = I - tau(i)*v*v'
			avec v(1: i-1) = 0
			v(i) = 1
			v(i+1 : m) = A(i+1: m,i)
			le vecteur tau etant connu.
	input : les dimensions m et n de la matrice A 
			le tableau A contenant les elements de la matrice apres facto QR (donnee par dgeqrf)
			le tableau tau (donne par dgeqrf)
	output : matrice Q de dimension mxm, resultat de la factorisation QR
*/

double * extraireQ( int m, int n, double *A, double *tau ){

	int i,j,k,l;
	double *Q, *H1,*H2;
	double *v;

	Q = calloc(m*m,sizeof(*A));
	v = calloc(m,sizeof(*A));
	H1 = calloc(m*m,sizeof(*A));
	H2 = calloc(m*m,sizeof(*A));

	for (i=0; i<m;i++){
		for (j=0;j<m;j++){
			Q[i+j*m]=0.0;
			H1[i+j*m]=0.0;
			H2[i+j*m]=0.0;
		}
	}  

	// construction du premier vecteur v
	k=0;
	v[0]=1;
	for (i=1;i<m;i++){
		v[i] = A[i+m*k];
	}

	// construction du premier reflecteur H1
	for (i=0; i<m;i++){
		for (j=0;j<m;j++){
			H1[i+m*j]= ((i==j) - tau[0]*v[i]*v[j]) ;
		}
	}
	
	// construction des autres reflecteurs
	for (k=1;k<n;k++){
		// construction des vecteurs v
		for (i=0;i<k;i++){
			v[i] = 0.0;
		}
		v[k] = 1;
		for (i=k+1;i<m;i++){
			v[i] = A[i+m*k];
		}
		for (i=0; i<m;i++){
			for (j=0;j<m;j++){
				H2[i+m*j]= ((i==j) - tau[k]*v[i]*v[j]) ;
			}
		}
		// et produit des reflecteurs avec les precedents
		for (i=0; i<m;i++){
			for (j=0;j<m;j++){
				for (l = 0; l<m;l++){
					Q[i+j*m] = H1[i+l*m]*H2[l+j*m] + Q[i+j*m];
				}
			}
		}
		for (i=0;i<m;i++){
			for (j=0;j<m;j++){
				H1[i+j*m] = Q[i+j*m];
				Q[i+j*m] = 0.0;	
			}
		}
	}
	
	free((void *) Q);
	free((void *) H2);
	free((void *) v);
	
	return H1;
}
