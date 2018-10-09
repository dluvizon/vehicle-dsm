#include <time.h>
#include <stdbool.h>
#include <locale.h>

#include "thog.h"

static const bool thog_print_settings = false;

/* FIXME: implement the return of this function as a pointer to a newly
 * allocated struct_thog */
struct_thog load_settings(const char *settings)
{
	struct_thog s;
	FILE *file;
	int ret;

	file = fopen(settings, "r");
	if (NULL == file) {
		printf("ERROR: (FIXME) FILE COULD NOT BE OPENED\n");
		return s;
	}

	char *old_locale = strdup(setlocale(LC_ALL, NULL));
	setlocale(LC_ALL, "C");

	ret = fscanf(file, "%d", &s.nh);
	if (!ret)
		goto err;
	ret = fscanf(file, "%d", &s.ncx);
	if (!ret)
		goto err;
	ret = fscanf(file, "%d", &s.ncy);
	if (!ret)
		goto err;
	ret = fscanf(file, "%d", &s.bpc);
	if (!ret)
		goto err;
	ret = fscanf(file, "%d", &s.norm);
	if (!ret)
		goto err;
	ret = fscanf(file, "%s", s.wnorm);
	if (!ret)
		goto err;
	ret = fscanf(file, "%lf", &s.rad);
	if (!ret)
		goto err;
	ret = fscanf(file, "%s", s.grad);
	if (!ret)
		goto err;
	ret = fscanf(file, "%s", s.hmetric);
	if (!ret)
		goto err;
	ret = fscanf(file, "%s", s.weight_function);
	if (!ret)
		goto err;
	ret = fscanf(file, "%d", &s.debug);
	if (!ret)
		goto err;

	if (thog_print_settings) {
		printf("load_settings:\n");
		printf("\ts.nh=%d\n", s.nh);
		printf("\ts.ncx=%d\n", s.ncx);
		printf("\ts.ncy=%d\n", s.ncy);
		printf("\ts.bpc=%d\n", s.bpc);
		printf("\ts.norm=%d\n", s.norm);
		printf("\ts.wnorm=%s\n", s.wnorm);
		printf("\ts.rad=%.6f\n", s.rad);
		printf("\ts.grad=%s\n", s.grad);
		printf("\ts.hmetric=%s\n", s.hmetric);
		printf("\ts.weight_function=%s\n", s.weight_function);
		printf("\ts.debug=%d\n", s.debug);
	}

	s.deformable_weights = FALSE;
	s.noc = s.ncx * s.ncy;
	s.nob = s.noc * s.bpc;
	fclose(file);

	setlocale(LC_ALL, old_locale);
	free(old_locale);
	return s;
err:
	setlocale(LC_ALL, old_locale);
	free(old_locale);
	printf("ERROR: (FIXME) Bad return value from fscanf\n");
	return s;
}

/**/
double classify(const unsigned char *image, int nrows, int ncols,
		int x, int y, int w, int h, struct svm_model* model,
		double *prob, struct_thog sthog)
{

   int margin = (int)(0.15*h + 0.5);

   if ((y - margin) > 0 && (h + margin) < nrows) {
      y = y - margin;
      h = h + 2*margin;
   }

   unsigned char *sample = (unsigned char *)malloc((w*h)*sizeof(unsigned char));

   /*Cropping the candidate text region from the original image: */
   int i, j, k, l;
   for (j = y, k = 0; j < (y+h); j++, k++) {
      for (i = x, l = 0; i < (x+w); i++, l++) {
         sample[w * k + l] = image[ncols * j + i];
      }
   }
   
   double *descriptor = thog (sample, h, w, sthog);

   /*Converting the descriptor to the SVM model: */
   svm_node* s = (struct svm_node *)malloc((sthog.nob+1)*sizeof(struct svm_node));

   for (i = 0; i < sthog.nob; i++) {
      s[i].index = i+1;
      s[i].value = descriptor[i];
   }
   s[i].index = -1;
   s[i].value = -1;

   double predict_label = svm_predict_probability (model, s, prob);
   free(s);
   free(descriptor);
   free(sample);

   return predict_label;  
}

/**/
double* thog (unsigned char *image, int nrows, int ncols, struct_thog sthog) {

   /*T-HOG settings: */
   int new_height = sthog.nh;
   int number_of_cells_x = sthog.ncx;
   int number_of_cells_y = sthog.ncy;
   int bins_per_cell = sthog.bpc;
   int image_normalization = sthog.norm;
   char *image_normalization_weight = sthog.wnorm;
   double image_normalization_weight_radius = sthog.rad;
   char *gradient_option = sthog.grad; 
   char *histogram_normalization_metric = sthog.hmetric; 
   char *weight_function = sthog.weight_function;
   int deformable_weights = sthog.deformable_weights;
   int debug = sthog.debug;
   /*End*/

   int i;

   int safe_margin = 1;

   double noise = 0.03;

   double black_level = 0.02; /*assumed black level of image*/

   int image_logscale = FALSE;

   int number_of_cells = number_of_cells_x * number_of_cells_y;

   int number_of_bins = number_of_cells * bins_per_cell; /*HOG bins*/

   /*Image resizing: if the new height is negative the region is not resized.*/
   unsigned char *resized = NULL;

   /*Getting the resized dimensions: */
   int rwidth;

   int rheight = new_height; 

   if (new_height > 0) {
      resized = applylanczos (image, ncols, nrows, new_height, 1, &rwidth);
   }
   else { printf("error: wrong image dimension\n"); exit(1); }

   /*number of pixels of the resized image.*/
   int n = rwidth * rheight;

   //double *dnorm = (double *)malloc(n * sizeof(double));

   //double *dtheta = (double *)malloc(n * sizeof(double));
   
   double *grey = (double *)malloc(n * sizeof(double));

   for (i = 0; i < n; i++) { grey[i] = (double)(resized[i]);}

   /*Creating a matrix to hold the cells histogram.*/
   double *cells_histogram = (double *)malloc(number_of_bins * sizeof(double)); 
   for (i = 0; i < number_of_bins; i++) {
      cells_histogram[i] = 0;
   }
   if ( (rwidth <= 2) || (rheight <= 2) ) {
      printf("too small dimension\n");
      return cells_histogram;
   }

   /*Computing weights {x_weight, y_weight} to normalize the resized image.*/
   int x_weight_rad = -1, y_weight_rad = -1;

   double *x_weight = NULL, *y_weight = NULL;

   if (strcmp(image_normalization_weight,"Gauss") == 0) {
      /*Gaussian weights*/
      x_weight_rad = choose_gaussian_weight_size (image_normalization_weight_radius * rheight);
      y_weight_rad = choose_gaussian_weight_size (1.0 * rheight/3.0);

      x_weight = (double *)malloc((2 * x_weight_rad + 1) * sizeof(double));
      y_weight = (double *)malloc((2 * y_weight_rad + 1) * sizeof(double));

      compute_gaussian_weights (x_weight, x_weight_rad);
      compute_gaussian_weights (y_weight, y_weight_rad);
   }
   else if (strcmp(image_normalization_weight,"Binomial") == 0) {
      /*Binomial weights*/
      x_weight_rad = y_weight_rad = (int) image_normalization_weight_radius;

      x_weight = (double *)malloc((2 * x_weight_rad + 1) * sizeof(double));
      y_weight = (double *)malloc((2 * y_weight_rad + 1) * sizeof(double));

      compute_binomial_weights (x_weight, x_weight_rad);
      compute_binomial_weights (y_weight, y_weight_rad);
   }
   else {
      printf("error: choose a valid weight image normalization\n");
      exit(1);
   }
 
   /*Convert to grey scale*/

   if (debug) {
      /*Writing image norm*/
      //image_functions.write_pgm (grey, rwidth, rheight, 0.0, 1.0, "grey");
   }

   /*Image normalization*/
   double *grel = NULL;
   if (image_normalization) {
      if ( (x_weight != NULL) &&  (y_weight != NULL) ) {
         grel = normalize_grey_image (grey, rwidth, rheight, x_weight, x_weight_rad, y_weight, y_weight_rad, noise);
      }
      else {
         exit(1);
      }
   }
   else {
      for (i = 0; i < n; i++) { grel[i] = grey[i]; }
   }

   /*convert to log scale*/
   if (image_logscale) {
      convert_to_log_scale(grel, rwidth, rheight, black_level);
   }
 
   //get_mag_theta (grel, dnorm, dtheta, rwidth, rheight, safe_margin, noise, gradient_option, debug);

   int x, y;
  
   /*Getting the baselines: */
   double xleft = safe_margin;
   double xright = rwidth - safe_margin;
   double ybot = rheight - safe_margin - 1;
   double ytop = safe_margin;

   double **mwtx = alloc_dmatrix (number_of_cells_x, rwidth); 
   for (x = safe_margin; x < (rwidth - safe_margin); x++) {
      int cx;
      for (cx = 0; cx < number_of_cells_x; cx++) {
         mwtx[x][cx] = cell_weight (weight_function, number_of_cells_x, cx, x, xright, xleft);
      }
   }   

   double **mwty = alloc_dmatrix (number_of_cells_y, rheight); 
   for (y = safe_margin; y < (rheight - safe_margin); y++) {
      int cy;
      for (cy = 0; cy < number_of_cells_y; cy++) {
         mwty[y][cy] = cell_weight (weight_function, number_of_cells_y, cy, y, ybot, ytop);
      }
   }  

   double eps = sqrt(2.0)*noise; /*Assumed deviation of noise in gradient*/

   double eps2 = eps * eps;

   int* vbin = (int *)malloc(2 * sizeof(int)); 

   double* factor = (double *)malloc(2 * sizeof(double));
         
   double *grad = (double *)malloc(2 * sizeof(double));;

   for (x = safe_margin; x < (rwidth - safe_margin); x++) {

      for (y = safe_margin; y < (rheight - safe_margin); y++) {

         int position = y * rwidth + x;
      
         /*Computing image gradients*/
         gradient_simple (grel, rwidth, rheight, x, y, grad);

         /*Computing the gradient norm but return zero if too small*/
         double d2 = grad[0]*grad[0] + grad[1]*grad[1];

         double dnorm = 0.0;
         if (d2 > eps2) {
             dnorm = sqrt(d2 - eps2);
         }

         /*Computing the gradient direction.*/
         double dtheta = atan2(grad[1], grad[0]);

         if (dtheta < 0) {
            //dtheta[position] += Math.PI;
            dtheta += 2*M_PI;
         }
         
         get_bin_pos (bins_per_cell, dtheta, vbin, factor);

         /*Computing the cells histogram and fuzzy weights: */
         int cx, cy;
 
         for (cx = 0; cx < number_of_cells_x; cx++) {

            for (cy = 0; cy < number_of_cells_y; cy++) {

               int c_pos = cy * number_of_cells_x + cx;

               int bin_pos1 = c_pos * bins_per_cell + vbin[0];

               int bin_pos2 = c_pos * bins_per_cell + vbin[1];

               cells_histogram[bin_pos1] += (dnorm * mwtx[x][cx] * mwty[y][cy]) * factor[0];

               cells_histogram[bin_pos2] += (dnorm * mwtx[x][cx] * mwty[y][cy]) * factor[1];
            }
         }
      } 
   }

   free(grad);
   free(vbin);
   free(factor);

   /*Normalize the histogram of each cell to unit L1 or L2 norm: */
   double sum = 0.0;

   int bin;

   /*Normalization sum: */
   for (bin = 0; bin < number_of_bins; bin++) {
      if (L1) { sum += cells_histogram[bin]; }
      else { sum += cells_histogram[bin] * cells_histogram[bin]; }
   }
   
   double cell_norm = 0.0;
   if (L1) { cell_norm = sum + 1.0 * number_of_bins; }
   else { cell_norm = sqrt(sum + 1.0 * number_of_bins); }

   /*Descriptor normalization: */
   for (bin = 0; bin < number_of_bins; bin++) {
      cells_histogram[bin] = (float)(cells_histogram[bin]/cell_norm);
   }

   disalloc_dmatrix (mwtx, rwidth);
   disalloc_dmatrix (mwty, rheight);
   free(x_weight);
   free(y_weight);
   free(grel);
   //free(dnorm);
   //free(dtheta);
   free(grey);
   free(resized);

   return cells_histogram;
}

/**/
void get_bin_pos (int bins_per_cell, double dtheta, int *bin, double *factor) {
                
   /*Computing the bin according the gradient direction.*/
   int full_circ = 0;

   double period = (full_circ ? 2 * M_PI : M_PI);

   double a = bins_per_cell*(dtheta/period+1);
   int ia = (int)floor(a);
   double frac = a - (double)ia;

   bin[0] = (ia + bins_per_cell) % bins_per_cell;
   bin[1] = (ia + 1) % bins_per_cell;
   factor[0] = (1-frac);
   factor[1] = frac;

   assert ( (bin[0] >= 0) && (bin[0] < bins_per_cell));
   assert ( (bin[1] >= 0) && (bin[1] < bins_per_cell));
}

/**/
void get_mag_theta (
   double *grel, double *dnorm, double *dtheta, 
   int rwidth, int rheight, 
   int safe_margin, double noise, char *gradient_option, 
   int debug
)
{
   int n = rwidth * rheight;

   double eps = sqrt(2.0)*noise; /*Assumed deviation of noise in gradient*/

   double eps2 = eps * eps;

   int x, y;

   for (x = safe_margin; x < (rwidth - safe_margin); x++) {

      for (y = safe_margin; y < (rheight - safe_margin); y++) {

         int position = y * rwidth + x;

         /*Computing image gradients*/
         double *grad = (double *)malloc(2 * sizeof(double));;

         if (strcmp(gradient_option,"Sobel") == 0) {
             gradient_sobel (grel, rwidth, rheight, x, y, grad);
         }
         else if (strcmp(gradient_option,"Simple") == 0) {
             gradient_simple (grel, rwidth, rheight, x, y, grad);
         }
         else  {
             printf("Error: choose a valid gradient option\n");
             exit(1);
         }

         /*Computing the gradient norm but return zero if too small*/
         double d2 = grad[0]*grad[0] + grad[1]*grad[1];

         if (d2 <= eps2) {
             dnorm[position] = 0.0;
         }
         else {
             dnorm[position] = sqrt(d2 - eps2);
         }

         /*Computing the gradient direction.*/
         dtheta[position] = atan2(grad[1], grad[0]);

         if (dtheta[position] < 0) {
            //dtheta[position] += Math.PI;
            dtheta[position] += 2*M_PI;
         }
         free(grad);
      }
   }
}

/*Sobel gradient: */
void gradient_sobel (double *image, int width, int height, int x, int y, double *grad) {

   int position = y * width + x;

   double vmo = image[position - 1];
   double vpo = image[position + 1];
   double vom = image[position - width];
   double vop = image[position + width];

   double vmm = image[position - 1 - width];
   double vmp = image[position - 1 + width];
   double vpm = image[position + 1 - width];
   double vpp = image[position + 1 + width];

   grad[0] = (vpm + 2*vpo + vpp - vmm - 2*vmo - vmp)/8.0;
   grad[1] = (vmm + 2*vom + vpm - vmp - 2*vop - vpp)/8.0;
}

/*Simple gradient: */
void gradient_simple (double *image, int width, int height, int x, int y, double *grad) {

   int position_x = y * width + x;
   int position_y = y * width + x;

   if (x == 0) {
      position_x = y * width + (x + 1);
   }
   if (x == (width-1)) {
      position_x = y * width + (x - 1);
   }

   if (y == 0) {
      position_y = (y + 1) * width + x;
   }

   if (y == (height - 1)) {
      position_y = (y - 1) * width + x;
   }

   double kxm = image[position_x - 1];
   double kxp = image[position_x + 1];
   double kym = image[position_y - width];
   double kyp = image[position_y + width];

   grad[0] = (kxp - kxm)/2;
   grad[1] = (kyp - kym)/2;
}

/* Computes the cell weight factor for one axis {z} (x or y). Given: 
 * number of cells {ncz}, cell index {cz}, pixel index {z}, all along 
 * that axis. The option selects the weight type. Adjusting the 
 * interval by {r_bot, r_top}.*/
double cell_weight (char *weight_function, int ncz, int cz, int z, double zmax, double zmin) {

   if (ncz == 1) { return 1; }

   double z_star = (z - zmin)/(zmax - zmin);

   if (strcmp(weight_function,"Step") == 0) {
      return StepFunc (ncz, cz, z_star);
   }
   else if (strcmp(weight_function,"Bernstein") == 0) {
      return Bernstein (ncz - 1, cz, z_star);
   }
   else if (strcmp(weight_function,"Core") == 0) {
      return EdgeCore (ncz, cz, z_star);
   }
   else if (strcmp(weight_function,"Exp") == 0){
      return Exp (ncz, cz, z_star);
   }
   else {
      assert(0);
      exit(1);
      return 0;
   }
}

/*Divides the interval [0-1] into {n} equal parts and returns 
* 1.0 if {z} is in part number {k} (0..n-1), 0 otherwise. */
double StepFunc (int n, int k, double z) {
   assert((k >= 0) && (k < n));
   if (z <= 0) { return (k == 0 ? 1 : 0); }
   else if (z >= 1) { return (k == (n - 1) ? 1 : 0); }
   else {
      return ( (k <= z*n) && (z*n < k+1) ? 1 : 0);
   }
}

/*Computes the Bernstein polynomial of degree {n} and index {k} for the argument {z}.*/
double Bernstein (int n, int k, double z) {
   assert((k >= 0) && (k <= n));
   if (z <= 0) { return (k == 0 ? 1 : 0); }
   else if (z >= 1) { return (k == n ? 1 : 0); }
   else {
      double zmax = ((double)k)/((double)n);
      return BernsteinPoly(n,k,z)/BernsteinPoly(n,k,zmax);
   }
}

/**/
double BernsteinPoly (int n, int k, double z) {
   assert((k >= 0) && (k <= n));
   double res = 1.0;
   int i;
   for (i = 0; i < k; i++) {
      res = (res * (n - i))/(i+1)*z;
   }
   return res*pow(1-z,n-k);
}

/*An edge-core weight function. If {n == 1} returns 1, if (n == 2) returns 
 *weight 1.0 near the edges, or 1.0 in the core region depending on {k}*/
double EdgeCore (int n, int k, double z) {
   assert(n == 2);
   assert((k >= 0) && (k < n));
   if ( (z <= 0) || (z >= 1) ) { return (k == 0 ? 1 : 0); }
   else {
      double v = 4 * z * (1 - z);
      v = v*v;
      return (k==0? 1 - v: v);
   }
}

/*An exponential weight function: */
double Exp (int n, int k, double z) {

   double mu = 0.01;
   double sigma = 0.5;
   assert((k >= 0) && (k < n));
   double avg = -mu + (1 + 2 * mu)/(double)(n - 1)*k;
   double dev = sigma/n;
   return gaussian (z, avg, dev);
}

/**/
double gaussian (double z, double mu, double sigma) {
   return exp(- ((z - mu)*(z - mu))/(2 * sigma * sigma));
}

/**/
int choose_gaussian_weight_size (double dev) {
   double rwt = 2*dev;
   return (int)(ceil(rwt));
}

/**/
void compute_binomial_weights (double *weight, int rwt) {
   int nwt = 2*rwt+1;
   weight[0] = 1;
   int i, j;
   for (i = 1; i < nwt; i++) {
      weight[i] = 0.0;
      for (j = i; j >=1; j--) {
         weight[j] = (weight[j] + weight[j-1])/2;
      }
      weight[0] /= 2;
   }
}

/**/
void compute_gaussian_weights (double *weight, int rwt) {
   double eps = 0.01;
   double a = -log(eps);
   int nwt = 2*rwt+1;

   int i; 
   for (i = 0; i < nwt; i++) {
      double z = (i - rwt)/((double)rwt);
      weight[i] = (1-z*z)*exp(-a*z*z);
   }
}

/**/
double* get_gray_image (unsigned char *image, int w, int h) {
   double *out = (double *)malloc(w * h * sizeof(double));
   int x, y;
   for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
         int p = y * w + x;
         double r = image[3 * y * w + 3 * x + 0]; /*red*/
         double g = image[3 * y * w + 3 * x + 1]; /*green*/
         double b = image[3 * y * w + 3 * x + 2]; /*blue*/
         r /= 255.0;
         g /= 255.0;
         b /= 255.0;
         out[p] = 0.299*r + 0.587*g + 0.114*b;
         //printf("%d = %f %f %f = %f\n", p, r, g, b, out[p]);
      }
   }
   return out;
}

/**/
void convert_to_log_scale (double *grey, int w, int h, double eps) {
   int position;
   int n = w * h;
   for (position = 0; position < n; position++) {
      grey[position] = (log(grey[position] + eps) - log(eps))/(log(1+eps)-log(eps));
   }
}

/*Get a normalized image from a grey level image.*/
double *normalize_grey_image (double *grey, int w, int h, double *x_weight, int x_rwt, double *y_weight, int y_rwt, double noise) {

   double* grel = (double *)malloc(w * h * sizeof(double));

   double AVG, DEV;

   int x, y;
   for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
         int position = y * w + x;
         AVG = get_grey_avg (grey, w, h, x, y, x_weight, x_rwt, y_weight, y_rwt);
         DEV = get_grey_dev (grey, w, h, x, y, x_weight, x_rwt, y_weight, y_rwt, AVG, noise);
         grel[position] = (grey[position] - AVG)/(3*DEV)+0.5;
         if (grel[position] < 0) { grel[position] = 0.0; }
         else if (grel[position] > 1) { grel[position] = 1.0; }
      }
   }
   return grel;
}

/*Get the averaged pixel value weighted by normalizing window: */
double get_grey_avg (double *grey, int w, int h, int x, int y, double *x_weight, int x_rwt, double *y_weight, int y_rwt) {
   double sum_vwt = 0.0, sum_wt = 0.0;
   int dx, dy;
   for (dy = -y_rwt; dy <= y_rwt; dy++) {
      int y1 = y + dy;
      if ( (y1 >= 0) && (y1 < h) ) {
         for (dx = -x_rwt; dx <= x_rwt; dx++) {
            int x1 = x + dx;
            if ( (x1 >= 0) && (x1 < w) ) {
               int position = y1 * w + x1;
               double v = grey[position];
               double wt = x_weight[x_rwt+dx]*y_weight[y_rwt+dy];
               sum_vwt += v * wt;
               sum_wt += wt;
            }
         }
      }
   }
   return sum_vwt/sum_wt;
}

/*Get the deviation of a pixel given a normalizing window: */
double get_grey_dev (double *grey, int w, int h, int x, int y, double *x_weight, int x_rwt, double *y_weight, int y_rwt, double AVG, double noise) {
   double sum_v2wt = 0.0, sum_wt = 0.0;
   int dx, dy;
   for (dy = -y_rwt; dy <= y_rwt; dy++) {
      int y1 = y + dy;
      if ( (y1 >= 0) && (y1 < h) ) {
         for (dx = -x_rwt; dx <= x_rwt; dx++) {
            int x1 = x + dx;
            if ( (x1 >= 0) && (x1 < w) ) {
               int position = y1 * w + x1;
               double v = grey[position]-AVG;
               double wt = x_weight[x_rwt+dx]*y_weight[y_rwt+dy];
               sum_v2wt += v * v * wt;
               sum_wt += wt;
            }   
         }
      }
   }
   return sqrt(sum_v2wt/sum_wt + noise*noise);
}


