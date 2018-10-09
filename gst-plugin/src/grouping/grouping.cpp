#include "grouping.h"

gset *make_sets(int n)
{
	int v;

	gset *s = (gset *) malloc(n * sizeof(gset));
	for (v = 0; v < n; ++v) {
		s[v].parent = v;
		s[v].rank = 0;
	}

	return s;
}

/**/
int get_number_of_different_sets (gset subsets[], int n) {
   int sum = 0;
   int v;
   int *count = (int*)malloc(n*sizeof(int));
   for (v = 0; v < n; v++) { count[v] = 0; }
   for (v = 0; v < n; v++) { count[subsets[v].parent] = 1; }
   for (v = 0; v < n; v++) { if (count[v]) { sum++; } }  
   free(count);
   return sum;
}

/**/
int get_max_set (gset subsets[], int n) {
   if (n == 0) { return 0; }
   int max = subsets[0].parent;
   int v;
   for (v = 1; v < n; v++) { 
      if (subsets[v].parent > max) {
         max = subsets[v].parent; 
      }
   }
   return max+1;
}

/**/
int find (gset subsets[], int i) {
   /*Find root and make root as parent of i (path compression): */
   if (subsets[i].parent != i) {
      subsets[i].parent = find (subsets, subsets[i].parent);
   }
   return subsets[i].parent;
}

/**/
void Union (gset subsets[], int x, int y) {
   int xroot = find (subsets, x);
   int yroot = find (subsets, y);
   /*Attach smaller rank tree under root of high rank tree (Union by Rank): */
   if (subsets[xroot].rank < subsets[yroot].rank) {
      subsets[xroot].parent = yroot;
   }
   else if (subsets[xroot].rank > subsets[yroot].rank) {
      subsets[yroot].parent = xroot;
   }
   /*If ranks are same, then make one as root and increment its rank by one: */
   else{
      subsets[yroot].parent = xroot;
      subsets[xroot].rank++; 
   }
}

region *create_region(void)
{
	region *r = (region *) malloc(sizeof(region));

	r->np = 0;
	vector_init(&(r->xp));
	vector_init(&(r->yp));

	r->box[0][0] = INTMAX;
	r->box[0][1] = INTMIN;
	r->box[1][0] = INTMAX;
	r->box[1][1] = INTMIN;

	return r;
}

void append_pixel(region *r, int x, int y)
{
	if (x < r->box[0][0]) {
		r->box[0][0] = x;
	}
	if (x > r->box[0][1]) {
		r->box[0][1] = x;
	}
	if (y < r->box[1][0]) {
		r->box[1][0] = y;
	}
	if (y > r->box[1][1]) {
		r->box[1][1] = y;
	}
	vector_add(&(r->xp), (void *) (long) x);
	vector_add(&(r->yp), (void *) (long) y);
	r->np++;
}

/* ------------------------------------------------------------
Computes an unnormalized eigenvector of the symmetric 2x2 
matrix [[XX,XY],[XY,YY]] corresponding to the eigenvalue {ev}. */
double* compute_axis (double ev, double XX, double XY, double YY) {
   double ax = ev - YY, ay = XY;
   double bx = XY, by = ev - XX;
   double d = ax*bx + ay*by;
   double *u = (double *)malloc(2*sizeof(double));
   if (d >= 0) { u[0] = ax+bx; u[1] = ay+by; }
   else { u[0] = ax-bx; u[1] = ay-by; }
   double um = 2*sqrt(ev)/hypot(u[0], u[1]);
   u[0] = u[0]*um; u[1] = u[1]*um;
   return u;
}

/*  */
double max (double a, double b) {
   return (a > b ? a : b);
}

/*  */
double min (double a, double b) {
   return (a < b ? a : b);
}

/**
 * (Re)computes the coordinates of the segment's center of mass
 * {center} from the pixel coordinates {xp,yp} and betas {mp}.
 */
void recompute_barycenter(region *r)
{
	double s_xb = 0;
	double s_yb = 0;
	double s_b = 0;
	int k;

	assert(r->np > 0);

	for (k = 0; k < r->np; k++) {
		double mass = 1.0;	/* this->mp.at(k); */
		double x = (double) ((size_t) vector_get(&(r->xp), k) + 0.5);
		double y = (double) ((size_t) vector_get(&(r->yp), k) + 0.5);

		s_xb += x * mass;
		s_yb += y * mass;
		s_b += mass;
	}
	/*assert (s_b == this->mass);*/
	r->center[0] = (s_b == 0 ? 0.5 : s_xb/s_b);
	r->center[1] = (s_b == 0 ? 0.5 : s_yb/s_b);
}

/* ----------------------------------------------------------------------
  (Re)computes the semiaxes {axes} of the segment's ellipse of inertia, from 
  the pixel coordinates {xp,yp}, their betas {mp}, and their barycenter 
  {center}. Treats each pixel as a square rather than a point mass. 
  !!! Should assume a Gaussian photo sampling kernel... !!! */
void recompute_axes (region *r) {
   assert (r->np > 0);
   double s_xxb = 0;
   double s_xyb = 0;
   double s_yyb = 0;
   double s_b = 0;
   double mpix = 1.0/12.0; /* X or Y moment of a square with unit mass and unit side. */
   int k;
   for (k = 0; k < r->np; k++) {
      double mass = 1.0;//this->mp.at(k);
      double x = (double) ((size_t)vector_get(&(r->xp), k) + 0.5 - r->center[0]);
      double y = (double) ((size_t)vector_get(&(r->yp), k) + 0.5 - r->center[1]);
      s_xxb += (x*x + mpix)*mass;
      s_xyb += x*y*mass;
      s_yyb += (y*y + mpix)*mass;
      s_b += mass;
   }
   /*assert(s_b == this->mass);*/
   assert(s_b > 0);
   double XX = s_xxb/s_b;
   double XY = s_xyb/s_b;
   double YY = s_yyb/s_b;
   if (XY == 0) {
      /* Axis-aligned ellipse or circle: */
      double r0 = 2*sqrt(max(XX,YY));
      double r1 = 2*sqrt(min(XX,YY));
      r->axes[0][0] = r0; r->axes[0][1] =  0;
      r->axes[1][0] =  0; r->axes[1][1] = r1;
   } else {
      /* Skew ellipse. */
      /* Find the eigenvalues {ev0,ev1} of the moment matrix {{XX, XY},{XY.YY}}: */
      double T = XX + YY;        /* Trace of matrix. */
      double D = XX*YY - XY*XY;  /* Determinant. */
      double h2 = T*T/4 - D;
      //assert(h2 >= 0);
      double h = sqrt(h2);
      double ev0 = T/2 + h;
      double ev1 = T/2 - h;
      //assert(ev0 >= ev1);
      //assert(ev1 >= 0);
      /* Compute the eigenvectors {(ux,uy),(vx,vy)}: */
      double *u = compute_axis(ev0, XX, XY, YY);
      r->axes[0][0] = u[0];
      r->axes[0][1] = u[1];
      free(u); 
      double *v = compute_axis(ev1, XX, XY, YY);
      r->axes[1][0] = v[0];
      r->axes[1][1] = v[1];
      free(v); 
   }
}

/* ---------------------------------------------
  Returns the lengths of the two axes of inertia 
  (longest first) as a vector of two doubles: */
double* radii (region *r) {
   double *vet = (double *)malloc(2*sizeof(double));
   vet[0] = hypot(r->axes[0][0], r->axes[0][1]);
   vet[1] = hypot(r->axes[1][0], r->axes[1][1]);
   return vet;
}

/*  */
int connected (region *si, region *sj, double t1, double t2, double t3) {

   double si_h = si->box[1][1] - si->box[1][0];
   double si_w = si->box[0][1] - si->box[0][0];
   double cx_i = si->box[0][0] + (si->box[0][1] - si->box[0][0])/2.0;
   double cy_i = si->box[1][0] + (si->box[1][1] - si->box[1][0])/2.0;

   double sj_h = sj->box[1][1] - sj->box[1][0];
   double sj_w = sj->box[0][1] - sj->box[0][0];
   double cx_j = sj->box[0][0] + (sj->box[0][1] - sj->box[0][0])/2.0;
   double cy_j = sj->box[1][0] + (sj->box[1][1] - sj->box[1][0])/2.0;

   double h = min (si_h, sj_h);
   double dx = abs (cx_i - cx_j) - (si_w + sj_w)/2.0;
   double dy = abs (cy_i - cy_j);

   double hmin = min (si_h, sj_h);
   double hmax = max (si_h, sj_h);

   //printf("cxi = %f, cyi = %f, cxj = %f, cyj = %f, sih = %f, sjh = %f, h = %f, dx : %f, dy : %f, t3: %f, conta 1: = %f, conta: 2 = %f\n",
   //cx_i, cy_i, cx_j, cy_j, si_h, sj_h, h, dx, dy, t3, hmax/hmin, t3*h);

   //if ( (abs(si_h - sj_h) < t1 * h) && (dx < t2 * h) && (dy < t3 * h) ) {
   if ( (hmax/hmin < 1.5) && (dy < t3 * h) && (dx < t2 * h) ) {
     // printf("conectei\n");
      return 1;
   }
   else {
     // printf("nao conectei\n");
      return 0;
   }
}

/* Assumes that {src} is a binary image. */
Vector find_connected_components4 (unsigned char *src, int nrows, int ncols, double rmin, double rmax)
{
	int pmin = 0;
	Vector cc;		/* List of conected components. */
	vector_init(&cc);
	unsigned char *image = (unsigned char *) copy_uchar(src, nrows, ncols);

	if ((ncols > 0) && (nrows > 0)) {
		/* Image is not empty. */
		/* Scan {img} for valid pixels: */
		int xv, yv;

		for (yv = 0; yv < nrows; yv++) {
			for (xv = 0; xv < ncols; xv++) {
				int v = image[yv * ncols + xv];

				if (v > pmin) {
					/**
					 * Collect the 4-connected component
					 * that contains {xv,yv}.
					 */
					region *r = create_region();

					/**
					 * Insert the starting pixel {x,y} in
					 * the stack:
					 */
					append_pixel(r, xv, yv);

					/* Clear it: */
					image[yv * ncols + xv] = 0;

					/**
					 * Now look for neighbors of
					 * stacked pixels:
					 *
					 * Pixels {seg.{x,y}[0..nz-1]} have
					 * had neighbors checked.
					 */
					int nz = 0;
					while (nz < r->np) {
						/**
						 * Pop next unchecked pixel
						 * from queue, get its col {xz}
						 * and row {yz}:
						 */
						int xz = (long) vector_get(&(r->xp), nz);
						int yz = (long) vector_get(&(r->yp), nz);
						nz++;

						/* Check its four neighbors: */
						int sx, sy;
						for (sx = -1; sx <= +1; sx++) {
							for (sy = -1; sy <= +1; sy++) {
								int xu = xz + sx;
								int yu = yz + sy;
								if ((xu >= 0) && (xu < ncols) && (yu >= 0) && (yu < nrows)) {
									int u = image[yu * ncols + xu];
									if (u > pmin){
										/* Add {xu,yu} to gset,clear it: */
										append_pixel(r, xu, yu);
										image[yu * ncols + xu] = 0;
									}
								}
							}
						}
					}
					recompute_barycenter(r);
					recompute_axes(r);
					double *rad = radii(r);

					/* Geometric constraints parameters: */
					/* Small and big image regions are discarded here:*/
					int wbox = r->box[0][1] - r->box[0][0];
					int hbox = r->box[1][1] - r->box[1][0];
					if ((wbox < rmin) || (hbox < rmin) || (wbox > rmax) || (hbox > rmax)) {
						vector_free(&(r->xp));
						vector_free(&(r->yp));
						free(r);
					} else {
						vector_add(&cc, r);
					}
					free(rad);
				}
			}
		}
	}

	free(image);
	return cc;
}

/* Assumes that {src} is a binary image. */
Vector find_connected_components_working(unsigned char *src,
		int nrows, int ncols, double rmin, double rmax,
		double fmin, int xmin, int ymin, int xmax, int ymax)
{
	int pmin = 0;
	Vector cc;		/* List of conected components. */
	vector_init(&cc);

	unsigned char *image = (unsigned char *) copy_uchar(src, nrows, ncols);

	if ((ncols > 0) && (nrows > 0)) {
		/* Image is not empty. */
		/* Scan {img} for valid pixels: */
		int xv, yv;

		for (yv = ymin; yv < ymax; yv++) {
			for (xv = xmin; xv < xmax; xv++) {
				int v = image[yv * ncols + xv];

				if (v > pmin) {
					/**
					 * Collect the 4-connected component
					 * that contains {xv,yv}.
					 */
					region *r = create_region();
					/** Insert the starting pixel {x,y}
					 * in the stack:
					 */
					append_pixel(r, xv, yv);
					/* Clear it: */
					image[yv * ncols + xv] = 0;

					/**
					 * Now look for neighbors of stacked
					 * pixels:
					 */
					int nz = 0;
					/** Pixels {seg.{x,y}[0..nz-1]} have
					 * had neighbors checked.
					 */
					while (nz < r->np) {
						/** Pop next unchecked pixel
						 * from queue, get its
						 * col {xz} and row {yz}:
						 */
						int xz = (int) (size_t) vector_get(&(r->xp), nz);
						int yz = (int) (size_t) vector_get(&(r->yp), nz);
						nz++;
						/* Check its four neighbors: */
						int sx, sy;
						for (sx = -1; sx <= +1; sx++) {
							for (sy = -1; sy <= +1; sy++) {
								int xu = xz + sx;
								int yu = yz + sy;

								if ((xu >= 0) && (xu < ncols) && (yu >= 0) && (yu < nrows)) {
									int u = image[yu * ncols + xu];
									if (u > pmin) {
										/* Add {xu,yu} to gset,clear it: */
										append_pixel(r, xu, yu);
										image[yu * ncols + xu] = 0;
									}
								}
							}
						}
					}
					recompute_barycenter(r);
					recompute_axes(r);
					double *rad = radii(r);

					/* Geometric constraints parameters: */
					/* Small and big image regions are discarded here: */
					/*
					double box_h = r->box[1][1] - r->box[1][0];
					double box_w = r->box[0][1] - r->box[0][0];
					*/
					if ((rad[0] >= rmin) && (rad[0] <= rmax) && (rad[1] >= fmin*rad[0])) {
						vector_add(&cc, r);
					} else {
						vector_free(&(r->xp));
						vector_free(&(r->yp));
						free(r);
					}
					free(rad);
				}
			}
		}
	}
	free(image);

	return cc;
}

Vector** make_chains(Vector cc, double t1, double t2, double t3, int *nsets)
{
	int i, j;

	gset *s = make_sets(vector_count(&cc));

	for (i = 0; i < vector_count(&cc); i++) {
		region *ri = (region *) vector_get(&cc, i);
		for (j = i+1; j < vector_count(&cc); j++) {
			region *rj = (region *) vector_get(&cc, j);
			if (connected(ri, rj, t1, t2, t3) && (find (s,i) != find (s,j))) {
				Union(s, find (s,i), find (s,j));
			}
		}
	}

	if (vector_count(&cc) == 0) {
		*nsets = 0;
	} else {
		*nsets = get_max_set(s, vector_count(&cc));
	}

	Vector **chains = (Vector **) malloc((*nsets) * sizeof(Vector *));

	for (i = 0; i < (*nsets); i++) {
		chains[i] = (Vector *) malloc(sizeof(Vector));
		vector_init(chains[i]);
	}

	for (i = 0; i < vector_count(&cc); i++) {
		region *ri = (region *) vector_get(&cc, i);
		vector_add(chains[find(s,i)], ri);
	}

	free(s);
	return chains;
}

Vector ** merge_chains  (Vector **chaina, int sizea, Vector **chainb, int sizeb) {
   int i, j;
   Vector **out = (Vector **)malloc((sizea + sizeb)*sizeof(Vector *));
   for (i = 0; i < sizea; i++) {
      out[i] = (Vector *)malloc(sizeof(Vector));
      vector_init(out[i]);
      for (j = 0; j < vector_count(chaina[i]); j++) {
          region *r = (region *)vector_get(chaina[i], j);
          vector_add (out[i], r);
      }
      vector_free(chaina[i]);
      free(chaina[i]);
   }
   for (i = 0; i < sizeb; i++) {
      out[i+sizea] = (Vector *)malloc(sizeof(Vector));
      vector_init(out[i+sizea]);
      for (j = 0; j < vector_count(chainb[i]); j++) {
          region *r = (region *)vector_get(chainb[i], j);
          vector_add (out[i+sizea], r);
      }
      vector_free(chainb[i]);
      free(chainb[i]);
   }
   free(chaina);
   free(chainb);
   return out;
}

Vector** make_chains_and_merge(Vector ncc, Vector icc,
		double t1, double t2, double t3, int *nsets)
{
	int i, j;

	for (i = 0; i < vector_count(&icc); i++) {
		vector_add(&ncc, (region *) vector_get(&icc, i));
	}

	gset *s = make_sets(vector_count(&ncc));

	for (i = 0; i < vector_count(&ncc); i++) {
		region *ri = (region *) vector_get(&ncc, i);
		for (j = i+1; j < vector_count(&ncc); j++) {
			region *rj = (region *) vector_get(&ncc, j);
			if (connected(ri, rj, t1, t2, t3) &&
					(find(s,i) != find (s,j))) {
				Union(s, find (s,i), find (s,j));
			}
		}
	}

	*nsets = get_max_set(s, vector_count(&ncc));

	Vector **chains = (Vector **) malloc((*nsets) * sizeof(Vector *));

	for (i = 0; i < (*nsets); i++) {
		chains[i] = (Vector *) malloc(sizeof(Vector));
		vector_init(chains[i]);
	}

	for (i = 0; i < vector_count(&ncc); i++) {
		region *ri = (region *) vector_get(&ncc, i);
		vector_add(chains[find(s,i)], ri);
	}

	vector_free(&ncc);
	free(s);

	return chains;
}

