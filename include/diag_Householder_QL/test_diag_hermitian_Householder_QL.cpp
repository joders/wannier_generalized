#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <complex>
#include <fstream>
#define cd complex <double>
#define I cd(0,1.0)


#ifndef diag_symmetric_Householder_QL_defined
#include "diag_hermitian_Householder_QL.cpp" 
#endif


int main()
{
	int n=3;

	cd **B=new cd*[n];
	for (int m=0; m<n; m++)
		B[m]=new cd[n];

	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
			B[i][j]=cos(sqrt(2)*10*double(i+2))+100*sin(sqrt(3)*10*double(j+0.3))+0.21*(j+0.4)*(i-9)*I;
	}
//make hermitian
	for (int i=0; i<n; i++)
	{
		for (int j=i; j<n; j++)
			B[j][i]=conj(B[i][j]);
	}

	for (int i=0; i<n; i++) B[i][i]=real(B[i][i]);
	

	double *eigenvalues=new double[n];

/*
for (int i=0; i<n; i++)
		{
		for (int j=0; j<n; j++)
		{
            printf("%f ",B[i][j]);
            }
		printf("\n");
	}
	*/
	
	
	bool with_eigen_vecs=1;
	diag_hermitian_Householder_QL(B,  eigenvalues, n,  with_eigen_vecs);


/*
printf("eigenstates:\n");
for (int i=0; i<n; i++)
		{
		for (int j=0; j<n; j++)
		{
            printf("%f ",B[i][j]);
            }
		printf("\n");
		}

printf("\n\neigenvalues:\n");
for (int i=0; i<n; i++)
		{
            printf("%f\n",eigenvalues[i]);
		}
*/


	return 0;
}






