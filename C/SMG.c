
/*
*	Sparse Matrix Generator implementation
*	Author: Ariel Livshits 2018
*/

#include "SMG.h"

#define LIGHT_CAP 0.3
#define HEAVY_CAP 0.7


void swap(int* elem1, int* elem2) {
	int tmp = *elem1;
	*elem1 = *elem2;
	*elem2 = tmp;
}

int bernoulli_sample(double p) {
	struct timespec tv;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tv);
	srand48(tv.tv_nsec);
	return drand48() < p;	
}

int randInt(int range[2]) {
	struct timespec tv;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tv);
	srand48(tv.tv_nsec);	
	return (int)(drand48()*(1 + range[1] - range[0])) + range[0];
}

void Knuth_random_select(int array[], int size, int numRand) {
	
	if (numRand == 0) return;
		
	int chosen = 0;
	double rem_size, rem_numRand;
	
	for (int i = 0; i < size; i++) {
		rem_size = size - i;
		rem_numRand = numRand - chosen;
		if (bernoulli_sample(rem_numRand/rem_size)) {
			array[i] = 1;
			chosen++;
		}
	}
}

void shuffle(int array[], int size) {
	struct timespec tv;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tv);
	srand48(tv.tv_nsec);

	int i, j;
	if (size > 1)
		for (i = size - 1; i > 0; i--) {
			j = (int)(drand48()*(i + 1));
			swap(array + j, array + i);
		}

}

int** Sparse_Matrix_Generator(int size, int num_light_flows, int num_heavy_flows, int* light_range, int* heavy_range) {

	if (!size || (!light_range && num_light_flows) || (!heavy_range && num_heavy_flows))
		return NULL;

	double light_sum = 0;
	double heavy_sum = 0;
	
	int** arr = (int**)malloc(size * sizeof(int*));
	int* rand_idx = (int*)malloc(size * sizeof(int));

	for (int i = 0; i < size; i++) {
		arr[i] = (int*)calloc(size, sizeof(int));
		rand_idx[i] = i;
	}

	for (int num = 0; num < num_light_flows; num++) {
		shuffle(rand_idx, size);
		for (int i = 0; i < size; i++) {
			arr[i][rand_idx[i]] += randInt(light_range);
			light_sum += arr[i][rand_idx[i]];
		}
	}
	
	for (int num = 0; num < num_heavy_flows; num++) {
		shuffle(rand_idx, size);
		for (int i = 0; i < size; i++) {
			arr[i][rand_idx[i]] += randInt(heavy_range);
			heavy_sum += arr[i][rand_idx[i]];
		}
	}
	printf("skew: %f\n", light_sum/(light_sum + heavy_sum));
	free(rand_idx);
	return arr;
}


void Free_Sparse_Matrix(int ** arr, int size)
{
	if (!arr) return;
	for (int i = 0; i < size; i++)
		free(arr[i]);

	free(arr);
}
