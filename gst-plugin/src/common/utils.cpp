#include "utils.h"

void int_array_set(int *array, int c, size_t nmemb)
{
	while (nmemb--) {
		*array++ = c;
	}
}

/* */
double dmax (double a, double b) {
  if (a > b) { return a; } else { return b; };
}

/* */
double dmin (double a, double b) {
  if (a < b) { return a; } else { return b; };
}

/* */
double** alloc_dmatrix (int ncols, int nrows) {
   int i;
   double **matrix = (double **)malloc(nrows * sizeof(double *));
   for (i = 0; i < nrows; i++) {
      matrix[i] = (double *)malloc(ncols * sizeof(double));
   }
   return matrix;
}

/* */
void disalloc_dmatrix (double **matrix, int nrows) {
   int i;
   for (i = 0; i < nrows; i++) {
      free(matrix[i]);
   }
   free(matrix);
}

/* */
float*** alloc_ndim_float_vect (int dim1, int dim2, int dim3) {
   int i, j;
   float*** array = (float ***)malloc(dim1*sizeof(float**));
   for (i = 0; i< dim1; i++) {
       array[i] = (float **) malloc(dim2*sizeof(float *));
       for (j = 0; j < dim2; j++) {
          array[i][j] = (float *)malloc(dim3*sizeof(float));
       }
   }
   return array;
}
