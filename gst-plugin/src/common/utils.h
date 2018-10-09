#ifndef __UTILS_H
#define __UTILS_H

#include <stdlib.h>

void int_array_set(int *array, int c, size_t nmemb);

inline int imax(int a, int b)
{
	return (a > b) ? a : b;
}

inline int imin(int a, int b)
{
	return (a < b) ? a : b;
}

double** alloc_dmatrix (int ncols, int nrows);

void disalloc_dmatrix (double **matrix, int nrows);


#endif /* __UTILS_H */
