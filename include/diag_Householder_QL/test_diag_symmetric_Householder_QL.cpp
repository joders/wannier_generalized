#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <complex>
#include <fstream>
#define cd complex <double>



#ifndef diag_symmetric_Householder_QL_defined
#include "diag_symmetric_Householder_QL.cpp" //contains QL diagonalization algo with implicit shifts for tridiagonal matrices
#endif

void Householder_symmetric_to_tridiag_old(double **a, int n, double d[], double e[]);


int main()
{
	int n=30;

	double **B=new double*[n];
	for (int m=0; m<n; m++)
		B[m]=new double[n];

	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
			B[i][j]=cos(sqrt(2)*10*double(i+2))+100*sin(sqrt(3)*10*double(j+0.3));
	}


	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
			B[i][j]=B[max(i,j)][min(i,j)];
	}

	double *eigenvalues=new double[n];


	/*
for (int i=0; i<n; i++)
		{
		for (int j=0; j<n; j++)
		{
            printf("%f ",B[i][j]);
            }
		printf("\n");
		}*/

	bool sort=1;
	bool with_eigen_vecs=1;
	printf("\n\nsstart");


	double **A=new double*[n];
	for (int m=0; m<n; m++)
		A[m]=new double[n];

	for(int run=0; run<10000; run++)
	{


		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
				A[i][j]=B[i][j];
		}	

		diag_symmetric_Householder_QL(A,  eigenvalues, n,  sort, with_eigen_vecs);

	}


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

	printf("\n\nstop");
	//int st;cin>>st;

	return 0;
}






