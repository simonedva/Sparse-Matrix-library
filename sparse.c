/**
 * @file sparse.c
 * @brief Library for sparse matrixes
 *
 * @author Simone De Vecchi (<simonedva@gmail.com>)
 * @version 1.1
 * @since 1.0
 */

#include "sparse.h"
#include <stdio.h>
#include "math.h"

/**
 * @brief Generate the sparse matrix of the matrix associated with pointer in passed
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in Pointer to the first element of the matrix
 * @param m Number of rows of the matrix
 * @param n	Number of columns of the matrix
 *
 * @return 0 if errors occurred
 */
int generateSparse(elem_t* out, const double* in, const int m, const int n) {

	//check
	if (out == NULL || in == NULL || m== 0 || n==0) {
		return 0;
	}

	//the first element is m,n,nnz for construction
	out->i = m;
	out->j = n;
	out->value = m*n;

	//construction of sparse matrix
	int nnz = 0; //when cycle finish, in nnz i have number of non-zeros
	for (int i = 0; i < m;i++) {
		for (int j= 0; j < n;j++) {

			if (*(in+j+i*n) >= INFVALUE) {
				(out+nnz+1)->i = i;
				(out+nnz+1)->j = j;
				(out+nnz+1)->value = *(in+j+i*n);
				nnz++;
			}

		}
	}

	//can't know if in has correct m,n but i know m*n is the worst case
	if (nnz > m*n) {
		return 0;
	}

	//if preallocated memory by caller isn't enough return 0
	if (nnz > (int) out->value) {
		return 0;
	}

	//setting the new nnz
	out->value = nnz;

	return 1;
}

/**
 * @brief Multiplies two sparse matrixes and stores the result in the sparse matrix pointed by out
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in1 Pointer to the first element of the first sparse matrix to multiply
 * @param in2 Pointer to the first element of the second sparse matrix to multiply
 *
 * @return 0 if errors occurred
 */
int multiplySparse(elem_t* out, const elem_t* in1, const elem_t* in2) {

	//checking if matrixes are compatible
	if (in1->j != in2->i) {
		return 0;
	}

	//very used variables
	int nout_new = 0; //number of nnz in the output matrix
	int i;
	int j;
	double value1;
	double value2;
	int found;

	//very used elements
	elem_t curr1;
	elem_t curr2;
	elem_t curr3;

	for (int k=0; k < (int) in1->value; k++) {

		curr1 = *(in1+k+1);
		i = curr1.i;
		j = curr1.j;
		value1 = curr1.value;

		for (int h = 0; h < (int) in2->value; h++) {

			curr2 = *(in2+h+1);

			if (curr2.i == j) { //if this row index has non-zero element

				value2 = (in2+h+1)->value;

				found = 0;
				for (int c = 0; c < nout_new; c++) {
					//searching if elem at index (i,c) already exists
					curr3 = *(out+c+1);

					if (curr3.i == i && curr3.j == curr2.j) {
						(out+c+1)->value += value1*value2;
						found = 1;
						break;
					}

				}

				if (!found) { //if doesn't exist such position then create a new one
					(out+nout_new+1)->i = i;
					(out+nout_new+1)->j = curr2.j;
					(out+nout_new+1)->value = value1*value2;
					nout_new++;
				}

			}
		}

	}


	//now in out i have nout_new elements
	out->value = nout_new;

	deleteZerosSparse(out);

	return 1;
}

/**
 * @brief Stores in the sparse matrix pointed by out the result of the addition of the two sparse matrixes pointed by in1 and in2
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in1 Pointer to the first element of the first sparse matrix to add
 * @param in2 Pointer to the first element of the second sparse matrix to add
 *
 * @return 0 if errors occurred
 */
int addSparse(elem_t* out, const elem_t* in1, const elem_t* in2) {

	//checking if matrixes are compatible
	if (in1->i != in2->i && in1->j != in2->j) {
		return 0;
	}

	//out almost surely contains each element of in1
	if(!copySparse(out,in1)) {
		return 0; //avoid errors
	}

	int nout = (int) out->value;
	int n2 = (int) in2->value;
	int nout_new = (int) in1->value;

	elem_t curr;
	int notfound;
	for(int k=0; k < n2;k++) {

		curr = *(in2+k+1);

		//trying to find if there is already an element in this (i,j)
		notfound = 1; //flag to know if i need to add a new elem_t
		for(int l=0;l < nout_new; l++) {
			if (curr.i == (out+l+1)->i && curr.j == (out+l+1)->j) {

				//updating value
				(out+l+1)->value += curr.value;

				notfound = 0;
				break; //interrupts inner cycle
			}
		}

		//otherwise adds another elem_t to out
		if (notfound) {
			(out+nout_new+1)->i = curr.i;
			(out+nout_new+1)->j = curr.j;
			(out+nout_new+1)->value = curr.value;
			nout_new++;
		}
	}

	out->value = nout_new;

	deleteZerosSparse(out);

	return 1;
}

/**
 * @brief Stores in the sparse matrix pointed by out the same element of the sparse matrix pointed by in (nnz could be different)
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in Pointer to the first element of sparse matrix to copy
 *
 * @return 0 if errors occurred
 */
int copySparse(elem_t* out, const elem_t* in) {

	//check
	if (out->value < in->value) {
		return 0;
	}

	//copying element per element
	for(int i = 0; i < (int) (in->value) + 1; i++) {
		(out+i+1)->i = (in+i+1)->i;
		(out+i+1)->j = (in+i+1)->j;
		(out+i+1)->value = (in+i+1)->value;
	}

	return 1;
}

/**
 * @brief Stores in the matrix pointed by out the full matrix of the sparse matrix pointed by in
 *
 * @param out Pointer to the first element of the result matrix
 * @param m Number of rows of the matrix pointed by out passed
 * @param n Number of columns of the matrix pointed by out passed
 * @param in Pointer to the first element of the sparse matrix to expand
 *
 * @return 0 if errors occurred
 */
int fullSparse(double* out, const int m, const int n, const elem_t* in) {

	//gettin information from first element
	int in_n = in->i;
	int in_m = in->j;
	int in_nnz = (int) in->value;

	//check
	if (m != in_m || n != in_n) {
		return 0;
	}

	//init matrix
	for (int k=0; k<m*n; k++) {
		*(out+k) = 0;
	}

	//setting non-zeros elements
	for (int k=0; k<in_nnz+1; k++) {
		//gettin current element
		elem_t curr = *(in+k+1);
		*(out + curr.i*in_n + curr.j) = curr.value;
	}

	return 1;
}

/**
 * @brief transpose the sparse matrix in his own location
 *
 * @param matrix Pointer to the first element of sparse the matrix
 *
 * @return 0 if errors occurred
 */
int transposeSparse(elem_t* matrix) {

	//check
	if (matrix == NULL) {
		return 0;
	}

	//gettin information
	int nnz = (int) matrix->value;

	//transpose
	int tmp;
	for (int k = 0; k < nnz+1;k++) { //also 0-th element because dimension changes
		tmp = (matrix+k)->i;
		(matrix+k)->i = (matrix+k)->j;
		(matrix+k)->j = tmp;
	}

	return 1;
}

/**
 * @brief Delete elements which are approximately zero from the sparse matrix
 *
 * @param out Pointer to the first element of the sparse matrix
 *
 * @return 0 if errors occurred
 */
int deleteZerosSparse(elem_t* out) {

	//deleting element if they are approximately 0
	for (int k = 0; k < (int) out->value; k++) {
		if (k >= (int) out->value) {
			break;
		}

		int pos = k+1;

		if (fabs((out+pos)->value) < INFVALUE) {

			//gettin last element
			int nnz = (int) out->value;

			//swap last element with the element to delete
			(out+pos)->i = (out+nnz)->i;
			(out+pos)->j = (out+nnz)->j;
			(out+pos)->value = (out+nnz)->value;

			//update nnz
			nnz--;
			out->value = nnz;

			k--; //if the swapped element is also 0. need to be check
		}
	}

	return 1;
}

/**
 * @brief Print to standard output the matrix
 *
 * @param matrix Pointer to the first element of the sparse matrix
 *
 * @return 0 if errors occurred
 */
int printSparse(const elem_t* matrix) {
	if (matrix == NULL) {
		return 0;
	}

	printf("Sparse matrix %dx%d:\n",matrix->i,matrix->j);
	for (int k= 0;k<(int)matrix->value; k++) {
		printf("(%d,%d) = %f\n",(matrix+k+1)->i,(matrix+k+1)->j,(matrix+k+1)->value);
	}

	return 1;
}
