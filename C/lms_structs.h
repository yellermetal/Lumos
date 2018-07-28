
/*
 *	lms_structs header
 *	Author: Ariel Livshits 2018
 */

#ifndef _LMS_STRUCTS_H_
#define _LMS_STRUCTS_H_

#include "Resources.h"
#include "string.h"
#include "math.h"
#include <stdio.h>

/* -----------------------------------------------------------------------------------
*	Function: Make_Matrix
*	Input: matrix dimensions, 2-D int array
*	Output: sparse matrix data structure pointer of type 'lms_mat_ptr'
*	Description: makes an efficient data structure for storing a sparse matrix
------------------------------------------------------------------------------------ */
lms_mat_ptr Make_Matrix(int dim, int_ptr arr[]);

/* -----------------------------------------------------------------------------------
*	Function: Copy_Matrix
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: pointer to a new matrix
*	Description: makes a deep copy of the input
------------------------------------------------------------------------------------ */
lms_mat_ptr Copy_Matrix(lms_mat_ptr mat_ptr); 

/* -----------------------------------------------------------------------------------
*	Function: Free_Matrix
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: N/A
*	Description: frees data structure and dynamically allocated memory
------------------------------------------------------------------------------------ */
void Free_Matrix(lms_mat_ptr mat_ptr);

/* -----------------------------------------------------------------------------------
*	Function: getNonzeroElemsIndices
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: a dynamically-allocated list of indices of nonzero elements
*	Description: returns a list of nonzero elements in O(nz) time complexity
------------------------------------------------------------------------------------ */
int_ptr getNonzeroElemsIndices(lms_mat_ptr mat_ptr);

/* -----------------------------------------------------------------------------------
*	Function: getNonzeroElemsValues
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: a dynamically-allocated list of values of nonzero elements
*	Description: returns a list of nonzero elements in O(nz) time complexity
------------------------------------------------------------------------------------ */
dub_ptr getNonzeroElemsValues(lms_mat_ptr mat_ptr);

/* -----------------------------------------------------------------------------------
*	Function: subtract_config
*	Input: pointer to matrix of type 'lms_mat_ptr', pointer to configuration of type 'lms_config_ptr'
*	Output: N/A
*	Description: subtracts coeff from matrix at coordinates defined in matching
------------------------------------------------------------------------------------ */
void subtract_config(lms_mat_ptr mat_ptr, lms_config_ptr config_ptr);

/* -----------------------------------------------------------------------------------
*	Function: Make_Config
*	Input: pointer to matrix of type 'lms_mat_ptr', a matching outputed by 
*		   Hopcroft Karp algorithm, size of matching, and a user-defined function
*		   for calculating desired coefficient metric from an array of chosen integers.
*
*	Output: a config data structure containing the matching array and calculated coefficient.
*
*	Description: calculates coefficient from chosen values by HopcroftKarp, creates 
*				 a data structure to represent the configuration and subtracting
				 the value of the coefficient from relevent indices in the matrix.
------------------------------------------------------------------------------------ */
lms_config_ptr Make_Config(lms_mat_ptr mat_ptr, int_ptr matching, int match_size, double (*coeff_func)(dub_ptr,int));

/* -----------------------------------------------------------------------------------
*	Function: Free_ConfigList
*	Input: pointer to a list of configurations
*	Output: N/A
*	Description: Frees allocated memory.
------------------------------------------------------------------------------------ */
void Free_ConfigList(lms_config_ptr config_list);

/* -----------------------------------------------------------------------------------
*	Function: threshold
*	Input: pointer to matrix of type 'lms_mat_ptr', integer threshold
*	Output: N/A
*	Description: sets all elements lower than sclr to be 0.
------------------------------------------------------------------------------------ */
void threshold(lms_mat_ptr mat_ptr, double sclr);

/* -----------------------------------------------------------------------------------
*	Function: ceiling
*	Input: pointer to matrix of type 'lms_mat_ptr', positive integer threshold
*	Output: N/A
*	Description: sets all elements higher than sclr to be sclr.
------------------------------------------------------------------------------------ */
void ceiling(lms_mat_ptr mat_ptr, double sclr);

/* -----------------------------------------------------------------------------------
*	Function: subtract
*	Input: pointer to matrix of type 'lms_mat_ptr', double value
*	Output: N/A
*	Description: Subtracts a value from all nonzero elements in the matrix
------------------------------------------------------------------------------------ */
void subtract(lms_mat_ptr mat_ptr, double sclr);

/* -----------------------------------------------------------------------------------
*	Function: pow
*	Input: pointer to matrix of type 'lms_mat_ptr', integer value
*	Output: N/A
*	Description: Raises all nonzero elements in matrix to a power of sclr
------------------------------------------------------------------------------------ */
void pow(lms_mat_ptr mat_ptr, int sclr);

/* -----------------------------------------------------------------------------------
*	Function: min
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: integer
*	Description: returns minimum nonzero value in the matrix. Returns 0 if there isn't one.
------------------------------------------------------------------------------------ */
double min(lms_mat_ptr mat_ptr);

/* -----------------------------------------------------------------------------------
*	Function: max
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: integer
*	Description: returns maximum value in the matrix. Returns 0 if there isn't one.
------------------------------------------------------------------------------------ */
double max(lms_mat_ptr mat_ptr);

/* -----------------------------------------------------------------------------------
*	Function: is_zero
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: boolean value
*	Description: check if matrix is the zero matrix
------------------------------------------------------------------------------------ */
bool is_zero(lms_mat_ptr mat_ptr);

/* -----------------------------------------------------------------------------------
*	Function: to_neg_int_mat
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: dynamically allocated integer matrix 
*	Description: translate lms_struct into a negetive integer matrix
*   Note: Need to know size of matrix to properly release this memory
------------------------------------------------------------------------------------ */
int_ptr* to_neg_int_mat(lms_mat_ptr mat_ptr);

/* -----------------------------------------------------------------------------------
*	Function: print_mat
*	Input: pointer to matrix of type 'lms_mat_ptr'
*	Output: N/A
*	Description: prints the matrix
------------------------------------------------------------------------------------ */
void print_mat(lms_mat_ptr mat_ptr);




#endif /*_LMS_STRUCTS_H_*/
