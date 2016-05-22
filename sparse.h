/**
 * @file sparse.h
 * @brief Library for sparse matrixes
 *
 * @author Simone De Vecchi (<simonedva@gmail.com>)
 * @version 1.1
 * @since 1.0
 */

#include <stdlib.h>
#include <stdio.h>

struct elem {
	int i;
	int j;
	double value;
};

typedef struct elem elem_t;

/* Numbers below this value are considered 0 and doesn't appear in the sparse matrix */
#define INFVALUE 0.001

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
int generateSparse(elem_t* out, const double* in, const int m, const int n);

/**
 * @brief Multiplies two sparse matrixes and stores the result in the sparse matrix pointed by out
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in1 Pointer to the first element of the first sparse matrix to multiply
 * @param in2 Pointer to the first element of the second sparse matrix to multiply
 *
 * @return 0 if errors occurred
 */
int multiplySparse(elem_t* out, const elem_t* in1, const elem_t* in2);

/**
 * @brief Multiplies a sparse matrix by a generic matrix and stores the result in the sparse matrix pointed by out
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in1 Pointer to the first element of the first sparse matrix to multiply
 * @param in2 Pointer to the first element of the second generic matrix to multiply
 * @param m Number of rows of the generic matrix
 * @param n Number of columns of the generic matrix
 *
 * @return 0 if errors occurred
 */
int multiplySparse_Matrix(elem_t* out, const elem_t* in1, const double* in2, const int m, const int n);

/**
 * @brief Stores in the sparse matrix pointed by out the result of the addition of the two sparse matrixes pointed by in1 and in2
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in1 Pointer to the first element of the first sparse matrix to add
 * @param in2 Pointer to the first element of the second sparse matrix to add
 *
 * @return 0 if errors occurred
 */
int addSparse(elem_t* out, const elem_t* in1, const elem_t* in2);

/**
 * @brief Stores in the sparse matrix pointed by out the same matrix pointed by in
 *
 * @param out Pointer to the first element of the result sparse matrix
 * @param in Pointer to the first element of sparse matrix to copy
 *
 * @return 0 if errors occurred
 */
int copySparse(elem_t* out, const elem_t* in);

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
int fullSparse(double* out, const int m, const int n, const elem_t* in);

/**
 * @brief transpose the sparse matrix in his own location
 *
 * @param matrix Pointer to the first element of sparse the matrix
 *
 * @return 0 if errors occurred
 */
int transposeSparse(elem_t* matrix);

/**
 * @brief Delete elements which are approximately zero from the sparse matrix
 *
 * @param out Pointer to the first element of the sparse matrix
 *
 * @return 0 if errors occurred
 */
int deleteZerosSparse(elem_t* out);

/**
 * @brief Print to standard output the matrix
 *
 * @param matrix Pointer to the first element of the sparse matrix
 *
 * @return 0 if errors occurred
 */
int printSparse(const elem_t* matrix);
