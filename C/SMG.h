

/*
*	Sparse Matrix Generator header
*	Author: Ariel Livshits 2018
*/

#ifndef _SMG_H_
#define _SMG_H_

/* Includes */
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int** Sparse_Matrix_Generator(int size, int light_matrices, int heavy_matrices, int* light_range, int* heavy_range);
void Free_Sparse_Matrix(int** arr, int size);

#endif /*_SMG_H_*/
