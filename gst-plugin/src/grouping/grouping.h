#ifndef __GROUPING_H
#define __GROUPING_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <queue>

#include "vector.h"
#include "image.h"

#define IMAX +9999999
#define IMIN -9999999

using namespace std;

typedef struct _pixel {
   int x;
   int y;
   int label;
} pixel;

typedef struct _region {
	int np;            /* Number of pixels in segment.                  */
	Vector xp;         /* List of {x} indices of pixels in segment.     */
	Vector yp;         /* List of {y} indices of pixels in segment.     */
	int box[2][2];     /* Columns and rows spanned by segment.          */
	double center[2];  /* Center of mass of pixels (weighted by {mass}. */
	double axes[2][2]; /* Semiaxes of the ellipse of inertia.           */
} region;

/* A structure to represent a subset for union-find: */
typedef struct _gset {
	int parent;
	int rank;
} gset;

gset* make_sets(int n);

int find(gset subsets[], int i);

void Union(gset subsets[], int x, int y);

region* create_region(void);

double* compute_axis(double ev, double XX, double XY, double YY);

double max(double a, double b);

double min(double a, double b);

void append_pixel(region *r, int x, int y);

int connected(region *si, region *sj, double t1, double t2, double t3);

void recompute_barycenter(region *r);

void recompute_axes(region *r);

double* radii(region *r);

int get_max_set(gset subsets[], int n);

Vector find_connected_components (unsigned char *src, unsigned char *msk, int nrows, int ncols, double rmin, double rmax, double fmin, int
xmin, int ymin, int xmax, int ymax);

Vector find_connected_components_working (unsigned char *src, int nrows, int ncols, double rmin, double rmax, double fmin, int xmin, int
ymin, int xmax, int ymax);

Vector find_connected_components3 (unsigned char *src, int nrows, int ncols, double rmin, double rmax);

Vector find_connected_components4 (unsigned char *src, int nrows, int ncols, double rmin, double rmax);

Vector** make_chains (Vector cc, double t1, double t2, double t3, int *nsets);

Vector ** merge_chains  (Vector **chaina, int sizea, Vector **chainb, int sizeb);

Vector find_connected_componentsB (
 unsigned char *src, int nrows, int ncols, double rmin, double rmax, double fmin, 
 int xmin, int ymin, int xmax, int ymax, 
 int stepx, int stepy );

Vector** make_chains_and_merge (Vector ncc, Vector icc, double t1, double t2, double t3, int *nsets);

#endif /* __GROUPING_H */
