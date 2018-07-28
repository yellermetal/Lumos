
/*
 *	lms_structs implementation
 *	Author: Ariel Livshits 2018
 */

#include "lms_structs.h"

#define END -1
#define NONZERO_ITERATOR(mat_ptr, iter) for( iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
void removeVal(lms_mat_ptr mat_ptr, int idx);

lms_mat_ptr Make_Matrix(int dim, int* arr[])
{
	if (arr == NULL) return NULL;

	lms_mat_ptr matrix = (lms_mat_ptr)malloc(sizeof(lms_mat_t));
	matrix->dim = dim;
	matrix->nz = 0;
	
	for (int row = 0; row < dim; row++)
		for (int col = 0; col < dim; col++)
			if (arr[row][col] > 0)
				matrix->nz++;

	matrix->nzo = matrix->nz;
	
	if (matrix->nz == 0) {
		matrix->elems = NULL;
		matrix->nz_arr = NULL;
		matrix->nz_head = NULL;
		return matrix;
	}

	matrix->elems = (mat_elem_ptr)malloc(matrix->nz * sizeof(mat_elem_t));
	matrix->nz_arr = (coord_ptr)malloc(matrix->nz * sizeof(coord_t));
	matrix->nz_head = matrix->nz_arr;

	int i = 0;
	for (int col = 0; col < dim; col++) {
		for (int row = 0; row < dim; row++) {
			if (arr[row][col] > 0) {
				matrix->elems[i].row = row;
				matrix->elems[i].col = col;
				matrix->elems[i].val = arr[row][col];
				matrix->nz_arr[i].prev = i - 1;
				matrix->nz_arr[i].next = i + 1;
				matrix->nz_arr[i].index = i;
				i++;
			}
		}
	}
	matrix->nz_arr[matrix->nz - 1].next = END;
	return matrix;
}

lms_mat_ptr Copy_Matrix(lms_mat_ptr mat_ptr)
{
	if (!mat_ptr) return NULL;

	lms_mat_ptr mat_cpy = (lms_mat_ptr)malloc(sizeof(lms_mat_t));
	memcpy(mat_cpy, mat_ptr, sizeof(lms_mat_t));

	mat_cpy->elems = (mat_elem_ptr)malloc(mat_ptr->nzo * sizeof(mat_elem_t));
	memcpy(mat_cpy->elems, mat_ptr->elems, mat_ptr->nzo * sizeof(mat_elem_t));

	mat_cpy->nz_arr = (coord_ptr)malloc(mat_ptr->nzo * sizeof(coord_t));
	memcpy(mat_cpy->nz_arr, mat_ptr->nz_arr, mat_ptr->nzo * sizeof(coord_t));

	int head_idx = mat_ptr->nz_head->index;
	mat_cpy->nz_head = mat_cpy->nz_arr + head_idx;

	return mat_cpy;
}

void Free_Matrix(lms_mat_ptr mat_ptr) 
{
	if (!mat_ptr) return;

	free(mat_ptr->elems);
	free(mat_ptr->nz_arr);
	free(mat_ptr);
	return;
}

int_ptr getNonzeroElemsIndices(lms_mat_ptr mat_ptr)
{
	if (!mat_ptr || (mat_ptr->nz ==0)) return NULL;
	
	int_ptr nz_elems = (int_ptr)malloc(sizeof(int) * mat_ptr->nz);
	int iter;
	int nzi = 0;
	for( iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		nz_elems[nzi++] = iter;

	return nz_elems;
}

dub_ptr getNonzeroElemsValues(lms_mat_ptr mat_ptr)
{
	if (!mat_ptr || (mat_ptr->nz ==0)) return NULL;
	
	dub_ptr nz_elems = (dub_ptr)malloc(sizeof(double) * mat_ptr->nz);
	int iter;
	int nzi = 0;
	for( iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		nz_elems[nzi++] = mat_ptr->elems[iter].val;

	return nz_elems;
}

lms_config_ptr Make_Config(lms_mat_ptr mat_ptr, int_ptr matching, int match_size, double (*coeff_func)(dub_ptr,int))
{
	if (!matching || !mat_ptr || !coeff_func) return NULL;

	dub_ptr match_values = (dub_ptr)malloc(match_size * sizeof(double));
	int iter, idx = 0;
	for( iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		if (mat_ptr->elems[iter].row == matching[mat_ptr->elems[iter].col]) 
			match_values[idx++] = mat_ptr->elems[iter].val;
	
	double coeff = coeff_func(match_values,match_size);
	
	free(match_values);

	lms_config_ptr config = (lms_config_ptr)malloc(sizeof(lms_config_t));
	config->coeff = coeff;
	config->matching = matching;
	config->next = NULL;
	
	subtract_config(mat_ptr, config);

	return config;

}

void Free_ConfigList(lms_config_ptr config_list)
{
	if (!config_list) return;
	
	lms_config_ptr iter;
	iter = config_list->next;
	while ( config_list != NULL ) {
		free(config_list->matching);
		free(config_list);
		config_list = iter;
		if (iter != NULL)
			iter = iter->next;
	}		
}

void threshold(lms_mat_ptr mat_ptr, double sclr)
{
	if (!mat_ptr) return;

	int iter;
	for( iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		if (mat_ptr->elems[iter].val < sclr) 
			removeVal(mat_ptr, iter);
}

void ceiling(lms_mat_ptr mat_ptr, double sclr)
{
	if (!mat_ptr || sclr == 0) return;

	int iter;
	for( iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		if (mat_ptr->elems[iter].val > sclr) 
			mat_ptr->elems[iter].val = sclr;
}

void subtract(lms_mat_ptr mat_ptr, double sclr)
{
	if (!mat_ptr) return;

	int iter;
	for (iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next) {
		if (mat_ptr->elems[iter].val > 0) {
			if (mat_ptr->elems[iter].val == sclr)
				removeVal(mat_ptr, iter);
			else
				mat_ptr->elems[iter].val -= sclr;
		}
	}
}

void subtract_config(lms_mat_ptr mat_ptr, lms_config_ptr config_ptr) 
{
	
	if (!mat_ptr || !config_ptr) return;
	
	double coeff = config_ptr->coeff;	
	int_ptr matching = config_ptr->matching;
	
	for(int iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		if (mat_ptr->elems[iter].row == matching[mat_ptr->elems[iter].col]) {
			mat_ptr->elems[iter].val -= coeff;
			if (mat_ptr->elems[iter].val <= 0)
				removeVal(mat_ptr, iter);
		}
	
}

void pow(lms_mat_ptr mat_ptr, int sclr)
{
	if (!mat_ptr) return;

	int iter;
	for (iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		mat_ptr->elems[iter].val = pow(mat_ptr->elems[iter].val, sclr);

}

double min(lms_mat_ptr mat_ptr)
{
	if (!mat_ptr) return 0;
	dub_ptr nz_elems = getNonzeroElemsValues(mat_ptr);
	double min = nz_elems[0];
	for (int i = 0; i < mat_ptr->nz; i++)
		if (nz_elems[i] < min)
			min = nz_elems[i];
	
	free(nz_elems);
	return min;
}

double max(lms_mat_ptr mat_ptr)
{
	if (!mat_ptr) return 0;

	int iter;
	double max = mat_ptr->elems[mat_ptr->nz_head->index].val;
	NONZERO_ITERATOR(mat_ptr, iter)
		if (mat_ptr->elems[iter].val > max)
			max = mat_ptr->elems[iter].val;

	return max;
}

bool is_zero(lms_mat_ptr mat_ptr)
{
	if (!mat_ptr) return true;
	return mat_ptr->nz == 0;
}

void removeVal(lms_mat_ptr mat_ptr, int idx) {

	int curr_next, curr_prev;
	mat_ptr->elems[idx].val = 0;
	mat_ptr->nz--;
	curr_next = mat_ptr->nz_arr[idx].next;
	curr_prev = mat_ptr->nz_arr[idx].prev;

	if (curr_next == END && curr_prev == END) {
		return;
	}

	if (curr_prev == END)
		mat_ptr->nz_head = mat_ptr->nz_arr + curr_next;
	else
		mat_ptr->nz_arr[curr_prev].next = curr_next;
	
	if (curr_next != END)
		mat_ptr->nz_arr[curr_next].prev = curr_prev;
	
}

int_ptr* to_neg_int_mat(lms_mat_ptr mat_ptr) {
	
	
	int size = mat_ptr->dim;
	
	int_ptr* matrix = (int_ptr*)malloc(size*sizeof(int_ptr));
	for (int i = 0; i < size; i++)
		matrix[i] = (int_ptr)malloc(size*sizeof(int));
	
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			matrix[i][j] = 0;
	
	for(int iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		matrix[mat_ptr->elems[iter].row][mat_ptr->elems[iter].col] = -1 * mat_ptr->elems[iter].val;

	return matrix;
}

void print_mat(lms_mat_ptr mat_ptr) {
	
	
	int size = mat_ptr->dim;
	
	dub_ptr* matrix = (dub_ptr*)malloc(size*sizeof(dub_ptr));
	for (int i = 0; i < size; i++)
		matrix[i] = (dub_ptr)malloc(size*sizeof(double));
	
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			matrix[i][j] = 0;
		
	for(int iter = mat_ptr->nz_head->index; iter != END; iter = mat_ptr->nz_arr[iter].next)
		matrix[mat_ptr->elems[iter].row][mat_ptr->elems[iter].col] = mat_ptr->elems[iter].val;
	
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			printf("%0.2f, ", matrix[i][j]);
		printf("\n");
		free(matrix[i]);
	}
	free(matrix);
	printf("\n");
	
}
