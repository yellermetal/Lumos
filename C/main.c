/*
 *Author: Ariel Livshits 2018
 */

/* Includes and Defines */
#include "lms_structs.h"
#include "HopcroftKarp.h"
#include "SMG.h"
#include "Algorithms.h"
#include <unistd.h>

#define TRACE_TYPE_SKEWED 1
#define TRACE_TYPE_STRONGLY_SKEWED 2
#define TRACE_TYPE_UNIFORM 3
#define TRACE_NUM 3
#define ALG_NUM 3
#define S2MS 1000
#define NS2MS 1000000

bool sanity_check(int** arr, lms_config_ptr config_list, int size);
double comp_time(lms_config_ptr config_list, double reconfiguration_delay);
void print_array(int** array, int size);
int** QuickStuff(int** matrix, int size);
bool check_DS(int** matrix, int size);
double avg(double array[], int size);
double variance(double array[], int size, double mean);
double gettime_ms(struct timespec start, struct timespec end);



int main(int argc, char *argv[]) {

	int ind = 1;
	int mat_num = atoi(argv[ind++]);
	int size = atoi(argv[ind++]);
	double delta = atoi(argv[ind++]);
	int light_matrices_num, heavy_matrices_num;
	int light_range[2], heavy_range[2];
	

	light_matrices_num = atoi(argv[ind++]);
	heavy_matrices_num = atoi(argv[ind++]);
	light_range[0] = atoi(argv[ind++]);
	light_range[1] = atoi(argv[ind++]);
	heavy_range[0] = atoi(argv[ind++]);
	heavy_range[1] = atoi(argv[ind++]);	

	if (mat_num < 0 || size < 0 || delta < 0) return 1;

	struct timespec t1[ALG_NUM], t2[ALG_NUM], QSt[2];
	lms_config_ptr  config, config_list[ALG_NUM];
	double elapsed_time, QS_elapsed_time;
	
	
	int** arr;
	int** sarr;
	int config_num;
	double* time[ALG_NUM];
	double* config_length[ALG_NUM];
	double* completion_time[ALG_NUM];
	const char* algo_names[ALG_NUM] = {"Solstice", "Eclipse", "Lumos"};
	
	for (int i = 0; i < ALG_NUM; i++) {
		time[i] = (double*)malloc(mat_num*sizeof(double));
		config_length[i] = (double*)malloc(mat_num*sizeof(double));
		completion_time[i] = (double*)malloc(mat_num*sizeof(double));
	}
	
	//char c;
	FILE* log = fopen("./logfile.txt", "w");
	FILE* time_log = fopen("./time_log.txt", "w");
	
	for (int i = 0; i < mat_num; i++) {
		
		arr = Sparse_Matrix_Generator(size, light_matrices_num, heavy_matrices_num, light_range, heavy_range);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &QSt[0]);
		sarr = QuickStuff(arr, size);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &QSt[1]);
		
		lms_mat_ptr mat_ptr = Make_Matrix(size, arr);
		lms_mat_ptr ds_mat_ptr = Make_Matrix(size, sarr);
		
		
		// Run and measure:
		int ind = 0;
		// SOLSTICE:
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1[ind]);
		config_list[ind] = SOLSTICE_Decomp(ds_mat_ptr);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2[ind++]);
		
		// ECLIPSE:
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1[ind]);
		config_list[ind] = ECLIPSE_Decomp(mat_ptr, delta);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2[ind++]);

		// LUMOS:
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1[ind]);
		config_list[ind] = LUMOS_Decomp(mat_ptr, delta);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2[ind++]);
		
		printf("cycle: %d, finished matrix num: %d\n", size, i + 1);

		
		for (int alg = 0; alg < ALG_NUM; alg++) {
			
			if(!sanity_check(arr, config_list[alg], size))
				printf("---%s failed.---\n", algo_names[alg]);
			
			// Count number of configurations:
			config_num = 0;
			for (config = config_list[alg]; config != NULL; config = config->next)
				config_num++;
			
			config_length[alg][i] = config_num;
			
			// Compute completion time:
			completion_time[alg][i] = comp_time(config_list[alg], delta);
		
			
			// Free heap memory:
			Free_ConfigList(config_list[alg]);
			

			// compute and print the elapsed time in millisec
			elapsed_time = gettime_ms(t1[alg], t2[alg]);   // ns to ms
			QS_elapsed_time = gettime_ms(QSt[0], QSt[1]);   // Quickstuff time
			

			time[alg][i] = elapsed_time;
			if (alg == 0) 
				time[alg][i] += QS_elapsed_time;
			
			fprintf(time_log, "% f, %f;", completion_time[alg][i], time[alg][i]);
			
		}
		
		fprintf(time_log, ";\n ");
		
		// Free heap memory:
		Free_Sparse_Matrix(arr, size);
		Free_Sparse_Matrix(sarr, size);
		Free_Matrix(mat_ptr);
		Free_Matrix(ds_mat_ptr);
		
	}

	double time_avg, time_var;
	double config_length_avg, config_length_var;
	double completion_time_avg, completion_time_var;
	
	for (int alg = 0; alg < ALG_NUM; alg++) {
		
		time_avg = avg(time[alg], mat_num);
		config_length_avg = avg(config_length[alg], mat_num);
		completion_time_avg = avg(completion_time[alg], mat_num);
		
		time_var = variance(time[alg], mat_num, time_avg);
		config_length_var = variance(config_length[alg], mat_num, config_length_avg);
		completion_time_var = variance(completion_time[alg], mat_num, completion_time_avg);


		fprintf(log, "%s %f %f %f %f %f %f \n" , algo_names[alg], config_length_avg, completion_time_avg, time_avg, config_length_var, completion_time_var, time_var);
		
		free(time[alg]);
		free(config_length[alg]);
		free(completion_time[alg]);
	}

	fclose(log);
	fclose(time_log);
	return 0;
}

double gettime_ms(struct timespec start, struct timespec end) {
	return (double)( (end.tv_sec - start.tv_sec)*S2MS + ((double)(end.tv_nsec - start.tv_nsec))/NS2MS );
}

double avg(double array[], int size) {
	
	double avg = 0;
	for (int i = 0; i < size; i++)
		avg += array[i];
	avg /= size;
	return avg;
}

double variance(double array[], int size, double mean) {
	double sqDiff = 0;
	for (int i = 0; i < size; i++)
		sqDiff += (array[i] - mean)*(array[i] - mean);
	sqDiff /= size;
	return sqDiff;
}

void print_array(int** array, int size) {

	for (int i = 0; i<size; i++) {
		for (int j = 0; j<size; j++)
			printf("%d, ", array[i][j]);
		printf("\n");
	}

	printf("\n");

}

double comp_time(lms_config_ptr config_list, double reconfiguration_delay) {
	
	int config_num = 0;
	int service_time = 0;
	lms_config_ptr iter = config_list->next;
	
	while (config_list != NULL) {
		
		service_time += config_list->coeff;
		config_num++;

		config_list = iter;
		if (iter != NULL)
			iter = iter->next;
	}
	return service_time + reconfiguration_delay * config_num;
}

bool sanity_check(int** arr, lms_config_ptr config_list, int size) {
	
	int** sarr = (int**)malloc(size * sizeof(int*));
	for (int i = 0; i < size; i++)
		sarr[i] = (int*)calloc(size, sizeof(int));

	lms_config_ptr iter;
	int coeff, *matching;

	iter = config_list->next;
	while (config_list != NULL) {

		matching = config_list->matching;
		coeff = config_list->coeff;

		for (int i = 0; i < size; i++) {
			if (matching[i] >= 0)
				sarr[matching[i]][i] += coeff;
		}

		config_list = iter;
		if (iter != NULL)
			iter = iter->next;
	}

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (arr[i][j] > sarr[i][j]) {				
				Free_Sparse_Matrix(sarr, size);
				return false;
			}
	
	Free_Sparse_Matrix(sarr, size);
	return true;
}

bool check_DS(int** matrix, int size) {
	
	int sum = 0, check_sum = 0;
	for (int j = 0; j < size; j++)
			sum += matrix[0][j];

	for (int i = 0; i < size; i++) {	
		check_sum = 0;
		for (int j = 0; j < size; j++)
			check_sum += matrix[i][j];
		
		if (check_sum != sum)
			return false;	
	}
	
	for (int j = 0; j < size; j++) {	
		check_sum = 0;
		for (int i = 0; i < size; i++)
			check_sum += matrix[i][j];
		
		if (check_sum != sum)
			return false;	
	}
	
	return true;
}

int** QuickStuff(int** matrix, int size) {
	
	int** stuffed_matrix = (int**)calloc(size, sizeof(int*));
	for(int i = 0; i < size; i++) {
		stuffed_matrix[i] = (int*)calloc(size, sizeof(int));
		memcpy(stuffed_matrix[i], matrix[i], size*sizeof(int));
	}
	
	int* row_sums = (int*)calloc(size, sizeof(int));
	int* col_sums = (int*)calloc(size, sizeof(int));
	
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			row_sums[i] += matrix[i][j];
			col_sums[j] += matrix[i][j];
		}
			
	int matrix_max = max(max(row_sums, size),max(col_sums, size));
	
	int stuffing = 0;
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (matrix[i][j] > 0) {
				stuffing = matrix_max - max(row_sums[i], col_sums[j]);
				stuffed_matrix[i][j] += stuffing;
				row_sums[i] += stuffing;
				col_sums[j] += stuffing;
			}
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (matrix[i][j] == 0) {
				stuffing = matrix_max - max(row_sums[i], col_sums[j]);
				stuffed_matrix[i][j] += stuffing;
				row_sums[i] += stuffing;
				col_sums[j] += stuffing;
			}	
	
	free(row_sums);
	free(col_sums);
	return stuffed_matrix;
}