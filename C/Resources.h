/*
*	Resources header
*	Author: Ariel Livshits 2018
*/

#ifndef _RESOURCES_H_
#define _RESOURCES_H_

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

typedef double* dub_ptr;
typedef int* int_ptr;
typedef bool* bool_ptr;

typedef struct _mat_elem_t {
	double val;
	int row;
	int col;
} mat_elem_t, *mat_elem_ptr;

typedef struct _coord_t {
	int next;
	int prev;
	int index;
} coord_t, *coord_ptr;

/* sols_mat_t is defines a sparse matrix structure */
typedef struct _lms_mat_t {

	mat_elem_ptr elems; /* all initial nonzero elements */
	coord_ptr nz_arr; /* indices of nonzero elements in elems */
	coord_ptr nz_head; /* pointer to first nonzero element in nz_arr */
	int nz; /* number of nonzero elements */
	int dim; /* the dimension of the matrix */
	int nzo; /* original nz */

} lms_mat_t, *lms_mat_ptr;

typedef struct _lms_config_t {
	double coeff; /* coefficient value */
	int_ptr matching; /* matching from HopcroftKarp */
	struct _lms_config_t * next;

} lms_config_t, *lms_config_ptr;

typedef struct _HopcroftKarp_t {
	int n, edges;
	int_ptr last, prev, head, matching, dist, Q;
	bool_ptr used, vis;
} HopcroftKarp_t, *HopcroftKarp_ptr;



#endif /*_RESOURCES_H_*/
