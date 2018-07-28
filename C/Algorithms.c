/*
*	BVN, Solstice, Eclipse & Lumos Algorithms header
*	Author: Ariel Livshits 2018
*/

/* Includes and Defines */

#include "Algorithms.h"

lms_config_ptr BVN_Decomp(lms_mat_ptr ds_mat_ptr) {


	lms_mat_ptr mat_ptr = Copy_Matrix(ds_mat_ptr);

	int_ptr matching = NULL;
	int size = 0;
	lms_config_ptr config = NULL, config_list = NULL;

	while (!is_zero(mat_ptr)) {

		HopcroftKarp(mat_ptr, &matching, &size);
		if (!config_list) {
			config_list = Make_Config(mat_ptr, matching, size, min);
			config = config_list;
		}
		else {
			config->next = Make_Config(mat_ptr, matching, size, min);
			config = config->next;
		}
	}
	Free_Matrix(mat_ptr);

	return config_list;

}

lms_config_ptr SOLSTICE_Decomp(lms_mat_ptr ds_mat_ptr) {


	lms_mat_ptr mat_ptr = Copy_Matrix(ds_mat_ptr);

	int_ptr matching = NULL;
	int size = 0;
	lms_config_ptr config = NULL, config_list = NULL;
	
	double th = pow(2,floor(log2(max(mat_ptr))));
	
	while (!is_zero(mat_ptr)) {
		
		lms_mat_ptr mat_ptr_temp = Copy_Matrix(mat_ptr);
		
		threshold(mat_ptr_temp, th);

		HopcroftKarp(mat_ptr_temp, &matching, &size);
		
		if (size == mat_ptr_temp->dim) {
			if (!config_list) {
				config_list = Make_Config(mat_ptr, matching, size, min);
				config = config_list;
			}
			else {
				config->next = Make_Config(mat_ptr, matching, size, min);
				config = config->next;
			}
		}
		else {
			free(matching);
			th /= 2;
		}
		matching = NULL;
		Free_Matrix(mat_ptr_temp);
	}

	Free_Matrix(mat_ptr);

	return config_list;

}

lms_config_ptr ECLIPSE_Decomp(lms_mat_ptr ds_mat_ptr, double delta) {


	lms_mat_ptr mat_ptr = Copy_Matrix(ds_mat_ptr);
	lms_config_ptr config = NULL, config_list = NULL;

	while (!is_zero(mat_ptr)) {
	
		if (!config_list) {
			config_list = Eclipse_Binary_Search(mat_ptr, delta);
			config = config_list;
		}
		else {
			config->next = Eclipse_Binary_Search(mat_ptr, delta);
			config = config->next;
		}
		
		subtract_config(mat_ptr, config);

	}
	Free_Matrix(mat_ptr);

	return config_list;

}

lms_config_ptr LUMOS_Decomp(lms_mat_ptr ds_mat_ptr, double delta) {

	lms_mat_ptr mat_ptr = Copy_Matrix(ds_mat_ptr);
	lms_config_ptr config = NULL, config_list = NULL;
		
	int mc = 0;
	int_ptr matching;

	while (!is_zero(mat_ptr)) {
	
		HopcroftKarp(mat_ptr, NULL, &mc);
		matching = Lumos_Binary_Search(mat_ptr, mc);
		
		
		if (!config_list) {
			config_list = Lumos_Config(mat_ptr, matching, mc, delta);
			config = config_list;
		}
		else {
			config->next = Lumos_Config(mat_ptr, matching, mc, delta);
			config = config->next;
		}

	}
	Free_Matrix(mat_ptr);

	return config_list;

}

lms_config_ptr Lumos_Config(lms_mat_ptr mat_ptr, int_ptr matching, int match_size, double delta) {
	
	if (!matching || !mat_ptr) return NULL;

	dub_ptr match_values = (dub_ptr)malloc(match_size * sizeof(double));
	int iter, idx = 0;
	for( iter = mat_ptr->nz_head->index; iter != -1; iter = mat_ptr->nz_arr[iter].next)
		if (mat_ptr->elems[iter].row == matching[mat_ptr->elems[iter].col]) 
			match_values[idx++] = mat_ptr->elems[iter].val;
	
	qsort(match_values, match_size, sizeof(double), dub_compare);
	
	double max_metric_value = -1;
	double normalized_edge_weight_sum;
	
	/*
	int start_idx = 0;
	for (int i = 0; i < match_size; i++)
		if (match_values[i] > 0) {
			start_idx = i;
			break;
		}
	*/
	double coeff = match_values[0];
	
	for (int i = 0; i < idx; i++) {
		
		normalized_edge_weight_sum = (i + 1) * (delta / (delta + match_values[i]));
		
		if (normalized_edge_weight_sum > max_metric_value) {
			max_metric_value = normalized_edge_weight_sum;
			coeff = match_values[i];
		}	
	}
	
	lms_config_ptr config = (lms_config_ptr)malloc(sizeof(lms_config_t));
	config->coeff = coeff;
	config->matching = matching;
	config->next = NULL;
	
	subtract_config(mat_ptr, config);
	free(match_values);
	return config;
}

int_ptr Lumos_Binary_Search(lms_mat_ptr ds_mat_ptr, int max_cardinality) {
	
	// get nonzero elems from matrix data-structure:
	dub_ptr elems_list = getNonzeroElemsValues(ds_mat_ptr);
	
	// unique sort the elements:
	int size = 0;
	
	if (ds_mat_ptr->nz == 0)
		return NULL;
	else
		elems_list = unique_sort(elems_list, ds_mat_ptr->nz, &size);
	
	if (size == 1) {
		double tmp = elems_list[0];
		free(elems_list);
		size = 2;
		elems_list = (dub_ptr)malloc(size*sizeof(double));
		elems_list[0] = tmp;
		elems_list[1] = tmp;
	}
	
	// Binary Search for :
	int mid , mc;
	int low = 0; 
	int high = size - 1;
	int_ptr matching;
	int_ptr last_good = NULL;
	
	lms_mat_ptr mat_ptr;
	
	while (low <= high) {
		
		mid = (low + high)/2;

		mat_ptr = Copy_Matrix(ds_mat_ptr);
		threshold(mat_ptr, elems_list[mid]);
		HopcroftKarp(mat_ptr, &matching, &mc);
		Free_Matrix(mat_ptr);

		if (mc == max_cardinality) {
			low = mid + 1;
			if (last_good != NULL) free(last_good);
			last_good = matching;
		}
		
		else {
			high = mid - 1;
			free(matching);
		}	
	}
	free(elems_list);	
	return last_good;
}

double Max_Weight_Matching(lms_mat_ptr ds_mat_ptr, double threshold, int_ptr matching) {
	
	if (!ds_mat_ptr || threshold == 0)	return 0;
	
	ceiling(ds_mat_ptr, threshold);	
	int_ptr* cost_matrix = to_neg_int_mat(ds_mat_ptr);
		
	hungarian_problem_t* hungarian_problem = (hungarian_problem_t*)malloc(sizeof(hungarian_problem_t));
	hungarian_init(hungarian_problem, cost_matrix, ds_mat_ptr->dim, ds_mat_ptr->dim, HUNGARIAN_MODE_MINIMIZE_COST);
	hungarian_solve(hungarian_problem);
	
	int_ptr* assignment = hungarian_problem->assignment;
	double Max_Weight = 0;
	
	for (int row=0; row < ds_mat_ptr->dim; row++)
		for (int col=0; col < ds_mat_ptr->dim; col++)
			if (assignment[row][col] == 1) {
				matching[col] = row;
				Max_Weight += -1 * cost_matrix[row][col];
			}
	
	Free_Matrix(ds_mat_ptr);
	hungarian_free(hungarian_problem);
	free(hungarian_problem);
	
	for (int i=0; i < ds_mat_ptr->dim; i++)
		free(cost_matrix[i]);
	free(cost_matrix);
	
	return Max_Weight;
	
}

lms_config_ptr Eclipse_Binary_Search(lms_mat_ptr ds_mat_ptr, double delta) {
	
	if (!ds_mat_ptr) return NULL;
	
	// get nonzero elems from matrix data-structure:
	dub_ptr elems_list = getNonzeroElemsValues(ds_mat_ptr);
	
	// unique sort the elements:
	int size = 0;
	
	if (ds_mat_ptr->nz == 0)
		return NULL;
	else
		elems_list = unique_sort(elems_list, ds_mat_ptr->nz, &size);
	
	if (size == 1) {
		double tmp = elems_list[0];
		free(elems_list);
		size = 2;
		elems_list = (dub_ptr)malloc(size*sizeof(double));
		elems_list[0] = tmp;
		elems_list[1] = tmp;
	}

	
	// Binary Search for Max Weight Matching using Hungarian Linear Assignment algorithm:
	int mid1, mid2;
	int low = 0; 
	int high = size - 1;
	double Max_Weight1 = 0;
	double Max_Weight2 = 0;
	int_ptr matching1 = (int_ptr)malloc(ds_mat_ptr->dim*sizeof(int));
	int_ptr matching2 = (int_ptr)malloc(ds_mat_ptr->dim*sizeof(int));
	
	lms_config_ptr config_ptr = (lms_config_ptr)malloc(sizeof(lms_config_t));
	config_ptr->next = NULL;
	
	while (low < high) {
		
		mid1 = (low + high)/2;
		mid2 = (low + high)/2 + 1;

		Max_Weight1 = Max_Weight_Matching(Copy_Matrix(ds_mat_ptr), elems_list[mid1], matching1) / (elems_list[mid1] + delta);
		Max_Weight2 = Max_Weight_Matching(Copy_Matrix(ds_mat_ptr), elems_list[mid2], matching2) / (elems_list[mid2] + delta);
		
		if (Max_Weight1 > Max_Weight2)
			high = mid1;
		
		else if( Max_Weight1 < Max_Weight2 )
			low = mid2;
		
		else {
			config_ptr->coeff =  elems_list[mid1];
			config_ptr->matching = matching1;
			free(matching2);
			free(elems_list);
			return config_ptr;
		}	
	}
	
	if (Max_Weight1 > Max_Weight2) {
		config_ptr->coeff =  elems_list[mid1];
		config_ptr->matching = matching1;
		free(matching2);
	}
	
	else {
		config_ptr->coeff =  elems_list[mid2];
		config_ptr->matching = matching2;
		free(matching1);
	}
	
	free(elems_list);
	return config_ptr;
	
}

int dub_compare (const void* a, const void* b) {
	return ( *(dub_ptr)a - *(dub_ptr)b );	
}

dub_ptr unique_sort(dub_ptr elem_list, int length, int_ptr new_length) {
	
	if (elem_list == NULL || length == 0)
		return 0;
	
	int size = length;
	qsort(elem_list, size, sizeof(double), dub_compare);
	
	for (int i = 1; i < length; i++)
		if (elem_list[i] == elem_list[i-1])
			size--;
	
	dub_ptr uqsort_elem_list = (dub_ptr)malloc(size*sizeof(double));
	
	int index = 1;
	uqsort_elem_list[0] = elem_list[0];
	for (int i = 1; i < length; i++)
		if (elem_list[i] != elem_list[i-1])
			uqsort_elem_list[index++] = elem_list[i];

	free(elem_list);
	*new_length = size;
	return uqsort_elem_list;	
}

void print_nz(lms_mat_ptr mat_ptr) {
	
	dub_ptr nz_elems = getNonzeroElemsValues(mat_ptr);
	for (int i = 0; i < mat_ptr->nz; i++)
		printf("%0.2f, ", nz_elems[i]);

	printf("\n");
	free(nz_elems);	
}

double min(dub_ptr elems, int size) {

	if (!elems || size == 0) return 0;

	double min = elems[0];
	for (int i = 0; i < size; i++)
		if (elems[i] < min)
			min = elems[i];

	return min;
}

double max(dub_ptr elems, int size) {

	if (!elems || size == 0) return 0;

	double max = elems[0];
	for (int i = 0; i < size; i++)
		if (elems[i] > max)
			max = elems[i];

	return max;
}

int min(int_ptr elems, int size) {

	if (!elems || size == 0) return 0;

	int min = elems[0];
	for (int i = 0; i < size; i++)
		if (elems[i] < min)
			min = elems[i];

	return min;
}

int max(int_ptr elems, int size) {

	if (!elems || size == 0) return 0;

	int max = elems[0];
	for (int i = 0; i < size; i++)
		if (elems[i] > max)
			max = elems[i];

	return max;
}

int min(int elem1, int elem2) {
	return (elem1 < elem2) ? elem1 : elem2;
}

int max(int elem1, int elem2) {
	return (elem1 > elem2) ? elem1 : elem2;
}

double mean(dub_ptr elems, int size) {

	if (!elems || (size == 0)) return 0;

	double avg = 0;
	for (int i = 0; i < size; i++)
		avg += elems[i];

	free(elems);
	return avg / size;
}

double Lumos_metric(dub_ptr elems, double val, int size) {
	
	double num_of_edges = 0, sum_of_negatives = 0;
	
	for (int i = 0; i < size; i++) {
		
		if (elems[i] <= val)	
			num_of_edges++;
		//else break;
		
		if (elems[i] < val)
			sum_of_negatives += (elems[i] - val);
		
	}
	
	return (num_of_edges + (sum_of_negatives/val))/size;
}

double Lumos_coeff_calc(dub_ptr elems, int size) {

	if (!elems || size == 0) return 0;

	//qsort(elems, size, sizeof(double), cmpfunc);
	
	double coeff = elems[0];
	for (int i = 0; i <size; i++)
		if ( Lumos_metric(elems, coeff, size) < Lumos_metric(elems, elems[i], size) )
			coeff = elems[i];
	

	return coeff;
}
