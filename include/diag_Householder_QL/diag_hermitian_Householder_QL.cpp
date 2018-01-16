#ifndef diag_hermitian_Householder_QL_defined
#define diag_hermitian_Householder_QL_defined

#ifndef diag_symmetric_Householder_QL_defined
#include "diag_symmetric_Householder_QL.cpp" //contains QL diagonalization algo with implicit shifts for tridiagonal matrices
#endif

void diag_hermitian_Householder_QL(cd **C,  double eigenvals[], int n,  bool with_eigen_vecs){
	// Transforms the problem of diagonalizing the n x n hermitian matrix into that of diagonalizing a 2n x 2n 
	// real, symmetric matrix using diag_symmetric_Householder_QL from the numerical recipes in C routine.
	// This may be a factor of 2 more inefficient than directly using a complex routine

	double eval_pair_tolerance=1E-10;

	//allocate 2n x 2n real matrix
	double **AB=new double*[2*n];
	for (int m=0; m<2*n; m++)	AB[m]=new double[2*n];

	cd **D=new cd*[n];	//temporary complex matrix
	for (int m=0; m<n; m++)	D[m]=new cd[2*n];
	
	cd *tmp_vec=new cd[n];
	
	for (int m1=0; m1<n; m1++){
		for (int m2=0; m2<n; m2++){
			AB[m1][m2]=real(C[m1][m2]);
			AB[n+m1][n+m2]=real(C[m1][m2]);
			AB[m1][n+m2]=-imag(C[m1][m2]);
			AB[n+m1][m2]=imag(C[m1][m2]);
		}
	}	
	
	/*
	
	FILE* F=fopen("input_matrix_tmp.m","w");
	fprintf(F,"\n\n\ninputmatrix=[\n");
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			//if(AB[i][j]>0) printf("+");
			fprintf(F,"%1.12f+(%1.12f*i) ",real(C[i][j]), imag(C[i][j]));
		}
		fprintf(F,"\n");
	}	
	fprintf(F,"];\n\n");
	fclose(F);
	*/
	double *AB_eigenvals=new double[2*n];
	int sort=1;
	diag_symmetric_Householder_QL(AB,  AB_eigenvals, 2*n,  sort, with_eigen_vecs);
	
	//check pairs
	for (int j=0; j<n; j++){	
		if(  fabs(AB_eigenvals[2*j]-AB_eigenvals[2*j+1])>eval_pair_tolerance  ) {
				printf("\nwarning: pair not found in subsequent eigenvals!!!!!!!!!!!!\n\n");
				exit(0);
		}
	}	
	
	int lo=0;
	int hi=0;
	int add_to_here;
	double norm;
	while(lo<2*n){
		hi=lo;
		while(fabs(AB_eigenvals[hi]-AB_eigenvals[lo])<eval_pair_tolerance && hi<2*n){hi++;}
		hi--;
		
		//now lo and hi are at the respective end of the same eigenvalue
		//check if degenerate
		if(  fabs(hi-lo)>1  ){//eigenvalue degenerate
			
			int num_vecs=(hi-lo+1)/2;
			//copy all complex vectors into temporary matrix
			for(int cmplx_vec_indx=0;cmplx_vec_indx<2*num_vecs;cmplx_vec_indx++){
				for(int pos=0; pos<n; pos++){D[pos][cmplx_vec_indx]=cd(AB[pos][lo+cmplx_vec_indx],AB[n+pos][lo+cmplx_vec_indx]);}
			}
			//now the overcomplete set of 2*num_vecs complex vectors have been saved to matrix D
			
			for(int j=0; j<num_vecs; j++){	//we are considering the jth vector to create linearly independent of the others
				
				add_to_here=j-1;
				do{
					add_to_here++;	//set to j for the first run
					for(int pos=0; pos<n; pos++) tmp_vec[pos]=D[pos][j];	//copy vector at position
					
					for(int j2=j+1;j2<=add_to_here;j2++){
						for(int pos=0; pos<n; pos++) tmp_vec[pos]+=D[pos][j2];
					}
					
					//subtract overlap of all previous vectors
					for(int j2=0; j2<j; j2++){ 
						cd overlap=0.0;
						for(int pos=0; pos<n; pos++) overlap+=conj(D[pos][j2])*tmp_vec[pos];		//calculate overlap 
						for(int pos=0; pos<n; pos++) tmp_vec[pos]-=overlap*D[pos][j2];	//subtract overlap*vector from tmp_vec
					}
					
					//normalize
					norm=0.0;
					for(int pos=0; pos<n; pos++) norm+=abs(tmp_vec[pos])*abs(tmp_vec[pos]);
					norm=sqrt(norm);
					//add so many later vectors until this is not the zero vector
				}while(norm<0.01 && add_to_here<(2*num_vecs-1));
				
				//now copy back to this column of D in normalized form
				for(int pos=0; pos<n; pos++) D[pos][j]=tmp_vec[pos]/norm;
			}
			
			//now the first num_vecs vectors in D still have to be copied to the correct position in the output C matrix
			for(int cmplx_vec_indx=0;cmplx_vec_indx<num_vecs;cmplx_vec_indx++){
				eigenvals[lo/2+cmplx_vec_indx]=AB_eigenvals[lo];	
				for (int pos=0; pos<n; pos++){
					C[pos][lo/2+cmplx_vec_indx]=D[pos][cmplx_vec_indx];
				}
			}
		}
		else{ //not degenerate
			eigenvals[lo/2]=AB_eigenvals[lo];
			for (int j2=0; j2<n; j2++){
				C[j2][lo/2]=cd( AB[j2][lo], AB[n+j2][lo]);
			}
		}
		//printf("lo:   EV[%d]=%f   hi:  EV[%d]=%f \n",lo,AB_eigenvals[lo],hi,AB_eigenvals[hi]);
		lo=hi+1;
	}
	
	/*
	for (int i=0; i<n; i++) printf("lambda[%d]=%1.14f\n",i,eigenvals[i]);
	printf("\n\n\noutput matrix with eigenvectors C:\n");
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			//if(AB[i][j]>0) printf("+");
			printf("%1.12f+(%1.12f*i) ",real(C[i][j]), imag(C[i][j]));
		}
		printf("\n");
	}	
	printf("\n\n");
	*/
	
	
	for (int m=0; m<n; m++) delete[] D[m];
	delete[] D;
	
	for (int m=0; m<2*n; m++) delete[] AB[m];
	delete[] AB;
	delete[] AB_eigenvals;
}

#endif //diag_hermitian_Householder_QL_defined
