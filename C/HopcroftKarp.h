
/*
*	HopcroftKarp header
*	Language: C
*	Author: Ariel Livshits 2018
*/

#ifndef _HOPCROFT_KARP_H_
#define _HOPCROFT_KARP_H_

/* Includes */
#include "Resources.h"
#include "lms_structs.h"

/*
 *	Function: HopcroftKarp algorithm
 *
 *	Input:    Matrix pointer of type 'lms_mat_ptr', an int pointer to hold matching array and a reference to int.
 *
 *	Output:   A dynamically-allocated array holding the matching found by HopcroftKarp algorthim.
 *			  The output array indices represent the columns and the values represent the matching rows.
 *			  The size of the array equals to the maximum cardinality matching and is outputed to out_maxcardinality.
 *
 *	Note:     If mat_ptr is NULL no output is generated. If max cardinality is zero, out_matching will be set to NULL,
			  and out_maxcardinality will be set to 0. Output "out_matching" must be passed to function calcCoefficient or freed.
 *
 *	Usage Example: 
 *		int size, *matching;
 *		lms_mat_ptr mat_ptr = Make_Matrix(dim, arr);
 *		HopcroftKarp(mat_ptr, matching, &size);
 *		int coeff = calcCoefficient(mat_ptr, matching, size, min);
 *		Free_Matrix(mat_ptr);
 *		
 */
void HopcroftKarp(lms_mat_ptr mat_ptr, int_ptr *out_matching, int_ptr out_maxcardinality );

#endif /*_HOPCROFT_KARP_H_*/