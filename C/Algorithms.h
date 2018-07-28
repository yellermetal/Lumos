/*
*	BVN, Solstice, Eclipse & Lumos Algorithms header
*	Author: Ariel Livshits 2018
*/

#ifndef _ALGORITHMS_H_
#define _ALGORITHMS_H_

#include "Resources.h"
#include "lms_structs.h"
#include "HopcroftKarp.h"
#include "SMG.h"
#include "hungarian.h"
#include <math.h>
#include <stdlib.h>

lms_config_ptr BVN_Decomp(lms_mat_ptr ds_mat_ptr);
lms_config_ptr SOLSTICE_Decomp(lms_mat_ptr ds_mat_ptr);
lms_config_ptr ECLIPSE_Decomp(lms_mat_ptr ds_mat_ptr, double delta);
lms_config_ptr LUMOS_Decomp(lms_mat_ptr ds_mat_ptr, double delta);
int_ptr Lumos_Binary_Search(lms_mat_ptr ds_mat_ptr, int max_cardinality);
lms_config_ptr Eclipse_Binary_Search(lms_mat_ptr ds_mat_ptr, double delta);
lms_config_ptr Lumos_Config(lms_mat_ptr mat_ptr, int_ptr matching, int match_size, double delta);
int dub_compare (const void* a, const void* b);
dub_ptr unique_sort(dub_ptr elem_list, int length, int_ptr new_length);
double Max_Weight_Matching(lms_mat_ptr ds_mat_ptr, double threshold, int_ptr matching);
double mean(dub_ptr elems, int size);
double Lumos_coeff_calc(dub_ptr elems, int size);
double min(dub_ptr elems, int size);
int min(int_ptr elems, int size);
double max(dub_ptr elems, int size);
int max(int_ptr elems, int size);
int min(int elem1, int elem2);
int max(int elem1, int elem2);

#endif /*_ALGORITHMS_H_*/
