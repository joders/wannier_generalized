#ifndef diag_symmetric_Householder_QL_defined
#define diag_symmetric_Householder_QL_defined

#include<cmath>
#include<cstdio>

static double sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
using namespace std;

//numerical recipes pg. 70
double pythag_ql(double a, double b) {
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);

	if (absa > absb)
		return absa * sqrt(1.0 + SQR(absb/absa));
	else
		return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa/absb)));
}

//**********************************************************************************************

void diag_tridiag_matrix_QL(double diagonal[], double offdiagonal[], int n, double *const*const z, bool sort, bool with_eigen_vecs)

/*
 in contrast to function in file diag_tridiag_matrix.cpp, this version does not transpose the output matrix z
 See numerical recipes in c, pg. 480
 QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric,
 tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2 §11.2. On
 input, diagonal[0..n-1] contains the diagonal elements of the tridiagonal matrix. On output, it returns
 the eigenvalues. The vector offdiagonal[0..n-2] inputs the subdiagonal elements of the tridiagonal matrix.
 On output offdiagonal is destroyed. When finding only the eigenvalues, several lines
 may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
 the matrix z[0..n-1][0..n-1] is input as the identity matrix. If the eigenvectors of a matrix
 that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
 In either case, the kth column of z returns the normalized eigenvector corresponding to diagonal[k].
 */
{
/*
printf("\n\ndiagonal elements in diag_tridiag_matrix_QL:   ");
for(int ctt=0; ctt<n; ctt++) printf("%f  ",diagonal[ctt]);

printf("\n\noffdiagonal elements in diag_tridiag_matrix_QL:   ");
for(int ctt=0; ctt<n-1; ctt++) printf("%f  ",offdiagonal[ctt]);
printf("\n\n\n");

exit(0);
*/


	int m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;

	for (l = 0; l < n; l++) {
		iter = 0;
		do {
			for (m = l; m < n - 1; m++) {
				dd = fabs(diagonal[m]) + fabs(diagonal[m + 1]);
				if (fabs(offdiagonal[m]) + dd == dd)
					break;
			}
			if (m != l) {
				if (iter++ == 40)
					printf("Too many iterations in diag_tridiag_matrix: iter=%d",iter);
				g = (diagonal[l + 1] - diagonal[l]) / (2.0 * offdiagonal[l]);
				r = pythag_ql(g, 1.0);
				g = diagonal[m] - diagonal[l] + offdiagonal[l]
						/ (g + SIGN(r,g));
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--) {
					f = s * offdiagonal[i];
					b = c * offdiagonal[i];
					offdiagonal[i + 1] = (r = pythag_ql(f, g));
					if (r == 0.0) {
						diagonal[i + 1] -= p;
						offdiagonal[m] = 0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = diagonal[i + 1] - p;
					r = (diagonal[i] - g) * s + 2.0 * c * b;
					diagonal[i + 1] = g + (p = s * r);
					g = c * r - b;
					// Next loop can be omitted if eigenvectors not wanted
					if (with_eigen_vecs) {
						for (k = 0; k < n; k++) {
							f = z[k][i + 1];
							z[k][i + 1] = s * z[k][i] + c * f;
							z[k][i] = c * z[k][i] - s * f;
						}
					}
				}
				if (r == 0.0 && i >= l)
					continue;
				diagonal[l] -= p;
				offdiagonal[l] = g;
				offdiagonal[m] = 0.0;
			}
		} while (m != l);
	}


	//sort eiganvals and eigenvecs if needed with heapsort according to code from numverical recipes in C
	// but indices from array always shifted fby 1 from original code, i.e. here the index of diagonal runs from 0 ... n-1
	if (sort) {
		int i, ir, j, l;
		double rra;
		double *zt = new double[n];

		if (n < 2)
			{
			delete[] zt;
			return;
			}
		l = (n >> 1) + 1;
		ir = n;
		for (;;) {
			if (l > 1) {
				rra = diagonal[--l-1];

				if (with_eigen_vecs) {
					for (int m = 0; m < n; m++)
						zt[m] = z[m][l-1];
				}

			} else {
				rra = diagonal[ir-1];

				if (with_eigen_vecs) {
					for (int m = 0; m < n; m++)
						zt[m] = z[m][ir-1];
				}

				diagonal[ir-1] = diagonal[0];

				if (with_eigen_vecs) {
					for (int m = 0; m < n; m++)
						z[m][ir-1] = z[m][0];
				}

				if (--ir == 1) {
					diagonal[0] = rra;

					if (with_eigen_vecs) {
						for (int m = 0; m < n; m++)
							z[m][0] = zt[m];
					}
					break;
				}
			}
			i = l;
			j = l + l;
			while (j <= ir) {
				if (j < ir && diagonal[j-1] < diagonal[j])
					j++;
				if (rra < diagonal[j-1]) {
					diagonal[i-1] = diagonal[j-1];

					if (with_eigen_vecs) {
						for (int m = 0; m < n; m++)
							z[m][i-1] = z[m][j-1];
					}

					i = j;
					j <<= 1;
				} else
					j = ir + 1;
			}
			diagonal[i-1] = rra;
			if (with_eigen_vecs) {
				for (int m = 0; m < n; m++)
					z[m][i-1] = zt[m];
			}
		}
		delete[] zt;
	}
}




void Householder_symmetric_to_tridiag(double *const*const a, int const n, double diagonal[], double offdiagonal[])
{
	/*
     From function "tred2.c" in numerical recipes in C. Householder transform for symmetric matrices, bringing these into
     tridiagonal form. Subsequently the tridiagonal matrix can be diagonalized by the QL algorithm with implicit shifts,
     defined in "diag_tridiag_matrix.cpp". However, the indices are shifted to C++ standard, i.e. A[0,...,n-1]


    Householder reduction of a real, symmetric matrix a[0..n-1][0..n-1]. On output, a is replaced
	by the orthogonal matrix Q effecting the transformation. diagonal[0..n-1] returns the diagonal elements
	of the tridiagonal matrix, and offdiagonal[0..n-1] the off-diagonal elements, with offdiagonal[0]=0. Several
	statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which
	case a contains no useful information on output. Otherwise they are to be included.
	 */
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n-1;i>=1;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				offdiagonal[i]=a[i][l];
			else {
				for (k=0;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				offdiagonal[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					offdiagonal[j]=g/h;
					f += offdiagonal[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<=l;j++) {
					f=a[i][j];
					offdiagonal[j]=g=offdiagonal[j]-hh*f;
					for (k=0;k<=j;k++)
						a[j][k] -= (f*offdiagonal[k]+g*a[i][k]);
				}
			}
		} else
			offdiagonal[i]=a[i][l];
		diagonal[i]=h;
	}
	diagonal[0]=0.0;

	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement diagonal[i]=a[i][i]; */
	for (i=0;i<n;i++) {
		l=i-1;
		if (diagonal[i]) {
			for (j=0;j<=l;j++) {
				g=0.0;
				for (k=0;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		diagonal[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}

	//shift to be compatible with convention used in diag_tridiag_matrix.cpp
	for (i=0;i<n-1;i++) offdiagonal[i]=offdiagonal[i+1];
	offdiagonal[n-1]=0.0;
}




/*************************************************************************************************
 **************************************************************************************************/




void diag_symmetric_Householder_QL(double *const*const a, double eigenvals[], int n,  bool sort, bool with_eigen_vecs)
{
	//according to numerical recipes, this is the fastest method for diagonalizing general real symmetric matrices
	//Input: symmetric matrix a[0..n-1][0..n-1]. On output, a is replaced by the unitary matrix of eigenvectors.
	//The eigenvalues are written into the array eigenvals[0..n-1]
	//the code is based on (but modified concerning the indices and incorporating sorting) the routes in "numerical recipes in C"

/*		FILE *F=fopen("Matrices.m","w");
	fprintf(F,"A=[\n");
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
			fprintf(F,"%1.14g ", a[i][j]);
		fprintf(F,"\n");
	}
	fprintf(F,"];\n\n");
*/

	double *offdiagonal=new double[n];
	Householder_symmetric_to_tridiag(a, n, eigenvals, offdiagonal);
	// a now contains the unitary matrix from the Householder transformation,
	// eigenvals the diagonal matrix elements and offdiagonal the off-diagonal elements of the tridiagonal matrix

	diag_tridiag_matrix_QL(eigenvals, offdiagonal, n, a, sort, with_eigen_vecs);
	delete[] offdiagonal;


	/*  //print out matrices in Matlab format for verification of algorithm
	fprintf(F,"V=[\n");
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
			fprintf(F,"%1.15g ", a[i][j]);
		fprintf(F,"\n");
	}
	fprintf(F,"];\n\n");

	fprintf(F,"\n\neigvals_cpp=[\n");
	for (int i=0; i<n; i++)
	{
		fprintf(F,"%f\n",eigenvals[i]);
	}
	fprintf(F,"]\n");
	fclose(F);
*/

}


#endif //diag_symmetric_Householder_QL_defined


