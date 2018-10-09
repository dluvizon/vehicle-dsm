/**
 * @file detection-zheng.cpp
 * @author Rodrigo Minetto
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 10/03/2015
 */

#define _DLEVEL SNOOPERTEXT_DLEVEL

#include <assert.h>
#include <limits.h>
#include <math.h>

#include "detection-zheng.h"
#include "grouping/grouping.h"
#include "seg/toggle.h"
#include "config.h"
#include "lanes.h"
#include "debug.h"

#define MAXGRAY 255

#define DEBUG 0

void write_plate (unsigned char *image, int ncols, int nrows, plate p, char *out_path, char *out_image_name, int region) {
   unsigned char *out = (unsigned char *)malloc(nrows * ncols * sizeof(unsigned char));
   memset (out, 0, nrows * ncols * sizeof(unsigned char));
   for (int y = p.y; y < (p.y + p.h); y++) {
      for (int x = p.x; x < (p.x + p.w); x++) {
         out[y * ncols + x] = image[y * ncols + x];
      }
   }
   write_uchar_to_pgm (out, ncols, nrows, "%s/%s_%02d_plate.pgm", out_path, out_image_name, region);
   free(out);
}

static double func_enhance (double std) {
   double r;
   if (std < 20) {
      r = 3.0 / (2.0/400.0 * (std - 20) * (std - 20) + 1);
   }
   else if (std < 60) {
      r = 3.0 / (2.0/1600.0 * (std - 20) * (std - 20) + 1);
   }
   else {
      r = 1.0;
   }
   return r;
} 

static void Enhancement (int xmin, int ymin, int xmax, int ymax, unsigned char *image, int nrows, int ncols, unsigned char *out, int wwin, int hwin) {

   xmin = max(xmin, wwin/2);
   ymin = max(ymin, hwin/2);
   xmax = min(xmax, ncols - wwin/2);
   ymax = min(ymax, nrows - hwin/2);

   for (int y = ymin; y <= ymax; y++) {
     for (int x = xmin; x <= xmax; x++) {
         double mean = 0.0;
         int sum = 0;
         for (int i = -hwin/2; i < +hwin/2; i++) {
            for (int j = -wwin/2; j < +wwin/2; j++) {
               int pixel = (y + i) * ncols + (x + j);
               mean += image[pixel];
               sum++;
            }
         }
         mean = mean/((double)(hwin*wwin));
         double std = 0.0;
         for (int i = -hwin/2; i < +hwin/2; i++) {
            for (int j = -wwin/2; j < +wwin/2; j++) {
               int pixel = (y + i) * ncols + (x + j);
               std += (image[pixel] - mean)*(image[pixel] - mean);
            }
         }
         std = sqrt(std/((double)(hwin*wwin)));
         int pixel = y * ncols + x;
         int value = (int)(func_enhance(std) * (image[pixel] - mean) + mean);
         if (value > 255) {
            out[pixel] = 255;
         }
         else if (value >= 0) {
            out[pixel] = (unsigned char)(value);
         }
         else {
            out[pixel] = image[pixel];
         }
      }
   }
}

/**/
void SobelZheng (int xmin, int ymin, int xmax, int ymax, unsigned char *image, int nrows, int ncols, unsigned char *out, int strenght) {

   int mask[3][3] = {{-1, +0, +1},
		     {-2, +0, +2},
		     {-1, +0, +1}};

   double *tmp = (double *)malloc(nrows * ncols * sizeof(double));
   
   int counter = 0;
 
   double sum = 0.0; 

   xmin = max(xmin, 1);
   ymin = max(ymin, 1);
   xmax = min(xmax, ncols - 1);
   ymax = min(ymax, nrows - 1);

   for (int y = ymin; y <= ymax; y++) {
      for (int x = xmin; x <= xmax; x++) {
         int i, j;
         double spixel = 0.0;
         for (j = -1; j <= 1; j++) {
	    for (i = -1; i <= 1; i++) {
               int pixel = (y + j)*ncols + (x + i);
	       spixel += mask[j+1][i+1] * image[pixel];
	    }
         }
         tmp[y * ncols + x] = fabs(spixel);
         sum += tmp[y * ncols + x];
         counter++;
      }
   }
   double threshold = strenght * (sum / double(counter));
   for (int y = ymin; y <= ymax; y++) {
      for (int x = xmin; x <= xmax; x++) {
         if (tmp[y * ncols + x] < threshold) {
            out[y * ncols + x] = 0;
         }
         else {
            out[y * ncols + x] = 255;
         }
      }
   }
   free(tmp); 
}

int maxA (int a, int b, int c, int d) {
   int v[4] = {a, b, c, d};
   int max = v[0];
   for (int i = 1; i < 4; i++) {
      if (v[i] > max) {
         max = v[i];
      }
   }
   return max;
}

int maxB (int a, int b, int c, int d, int e, int f) {
   int v[6] = {a, b, c, d, e, f};
   int max = v[0];
   for (int i = 1; i < 6; i++) {
      if (v[i] > max) {
         max = v[i];
      }
   }
   return max;
}

/**/
void Filtering (int xmin, int ymin, int xmax, int ymax, unsigned char *edges, int nrows, int ncols, int tshort, int tlong) {

   xmin = max(xmin, 2);
   ymin = max(ymin, 2);
   xmax = min(xmax, ncols - 2);
   ymax = min(ymax, nrows - 2);

   int *N = (int *)malloc(nrows * ncols * sizeof(int));
   
   int *M = (int *)malloc(nrows * ncols * sizeof(int));
   
   for (int x = 0; x < (nrows * ncols); x++) { 
      N[x] = M[x] = 0; 
   } 

   for (int i = ymin; i <= ymax; i++) {
      for (int j = xmin; j <= xmax; j++) {
         edges[i * ncols + j] /= 255;
      }
   }

   for (int i = ymin; i <= ymax; i++) {
      for (int j = xmin; j <= xmax; j++) {
         if (edges[i * ncols + j] == 1) {
            if ( (edges[(i-1) * ncols + (j-1)] + edges[(i-1) * ncols + (j)] + edges[(i-1) * ncols + (j+1)] + edges[(i) * ncols + (j-1)]) > 0) {
               M[i * ncols + j] = maxA (M[(i-1) * ncols + (j-1)], M[(i-1) * ncols + j], M[(i-1) * ncols + (j+1)], M[(i) * ncols + (j-1)]) + 1;
            }
            else {
               M[i * ncols + j] = maxB (M[(i-2) * ncols + (j-1)], M[(i-2) * ncols + j], M[(i-2) * ncols + (j+1)], 
                                     M[(i-1) * ncols + (j-2)], M[(i-1) * ncols + (j+2)], M[(i) * ncols + (j-2)]) + 1;
            }
         }
      }
   }
   for (int i = ymax; i >= ymin; i--) {
      for (int j = xmax; j >= xmin; j--) {
         if (edges[i * ncols + j] == 1) {
            if ( (edges[(i+1) * ncols + (j-1)] + edges[(i+1) * ncols + (j)] + edges[(i+1) * ncols + (j+1)] + edges[(i) * ncols + (j+1)]) > 0) {
               N[i * ncols + j] = maxA (N[(i+1) * ncols + (j-1)], N[(i+1) * ncols + j], N[(i+1) * ncols + (j+1)], N[(i) * ncols + (j+1)]) + 1;
            }
            else {
               N[i * ncols + j] = maxB (N[(i+2) * ncols + (j-1)], N[(i+2) * ncols + j], N[(i+2) * ncols + (j+1)], 
                                        N[(i+1) * ncols + (j-2)], N[(i+1) * ncols + (j+2)], N[(i) * ncols + (j+2)]) + 1;
            }
         }
      }
   }
   for (int i = ymin; i <= ymax; i++) {
      for (int j = xmin; j <= xmax; j++) {
         if (edges[i * ncols + j] == 1) {
            if ( (M[i * ncols + j] + N[i * ncols + j] > tlong) || (M[i * ncols + j] + N[i * ncols + j] < tshort) ) {
               edges[i * ncols + j] = 0;
            }
         }
      }
   }
   for (int i = ymin; i <= ymax; i++) {
      for (int j = xmin; j <= xmax; j++) {
         edges[i * ncols + j] *= 255;
      }
   }
   free(N); 
   free(M); 
}

/**/
plate find_plate (int xmin, int ymin, int xmax, int ymax, unsigned char *edge, int nrows, int ncols, int window_w, int window_h, double percentage)
{
   int *acumulator = (int *) malloc(nrows * ncols * sizeof(int));
   
   for (int i = 0; i < (nrows * ncols); i++) {
      acumulator[i] = 0;
   }

   xmin = max(xmin, window_w/2);
   ymin = max(ymin, window_h/2);
   xmax = min(xmax, ncols - window_w/2);
   ymax = min(ymax, nrows - window_h/2);

   for (int i = ymin; i <= ymax; i++) {
      for (int j = xmin; j <= xmax; j++) {
         int sum = 0; 
         for (int y = (i - window_h/2); y < (i + window_h/2); y++) {
            for (int x = (j - window_w/2); x < (j + window_w/2); x++) {
               if (edge[y * ncols + x] == 255) {
                  sum++; 
               }
            }
         }
         acumulator[i * ncols + j] = sum;
      }
   }

   int vmax = 0, xcenter = 0, ycenter = 0;
   for (int i = ymin; i <= ymax; i++) {
      for (int j = xmin; j <= xmax; j++) {
         if (acumulator[i * ncols + j] > vmax) {
            vmax = acumulator[i * ncols + j];
            xcenter = j; 
            ycenter = i; 
         }
      }
   }

   /*If the window has {percentage} of edge pixels then accept the region as a license plate.*/
   int threshold = (int)(window_w * window_h * percentage);
   if (vmax > threshold) {
      plate p = {xcenter - window_w/2, ycenter - window_h/2, window_w, window_h};
      p.speed = 0.0;
      p.from = NULL;
      free(acumulator);
      return p;
   }
   else {
      plate p = {0, 0, 0, 0};
      free(acumulator);
      return p;
   } 
}

vector<plate> detection_zheng (
	unsigned char *image,
	struct snoopertext *s,
	int ncols, int nrows,
	char *out_path,
	char *out_image_name,
	vector<struct sub_slope> &slopes,
	int margin_X, int margin_Y,
	int iframe,
	unsigned char *umask)
{
   clock_t start = clock();
   vector<plate> plates;

   /*Parameter settings:*/
      /*Enhancement: */
      int hwin = 8;             /*fixed*/
      int wwin = 8;             /*fixed*/
      /*Sobel: */
      int strenght = 3;         /*change*/
      /*Egde filtering: */
      int tshort = 16;           /*change*/
      int tlong = 68;           /*change*/
      /*Region filtering: */
      int window_w = 107;       /*change or FIXED?*/
      int window_h = 30;        /*change or FIXED?*/
      double percentage = 0.2;  /*change*/
   /*End settings*/

   int cter = 0; 
   for (int i = 0; i < slopes.size(); i++) {
      if (!slopes.at(i).search_plate) { continue; }
      cter++;
   }
   if (0 == cter) {
      return plates;
   }

   memset(s->edges, 0, nrows * ncols * sizeof(unsigned char));
   memset(s->enhanced, 0, nrows * ncols * sizeof(unsigned char));

   for (int i = 0; i < slopes.size(); i++) {
      if (!slopes.at(i).search_plate) {
         continue;
      }
      if (DEBUG) {
         write_uchar_to_pgm (image, ncols, nrows, "%s/%s_%02d_orig.pgm", out_path, out_image_name, i);
      }
      Enhancement (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, image, nrows, ncols, s->enhanced, wwin, hwin);
      if (DEBUG) {
         write_uchar_to_pgm (s->enhanced, ncols, nrows, "%s/%s_%02d_edges.pgm", out_path, out_image_name, i);
      }
      SobelZheng (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, s->enhanced, nrows, ncols, s->edges, strenght);
      if (DEBUG) {
         write_uchar_to_pgm (s->edges, ncols, nrows, "%s/%s_%02d_sobel.pgm", out_path, out_image_name, i);
      }
      Filtering (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, s->edges, nrows, ncols, tshort, tlong);
      if (DEBUG) {
         write_uchar_to_pgm (s->edges, ncols, nrows, "%s/%s_%02d_filtered.pgm", out_path, out_image_name, i);
      }
      plate p = find_plate (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, s->edges, nrows, ncols, window_w, window_h, percentage);
      if ( (p.w != 0) && (p.h != 0) ) {
         plates.push_back(p);
         if (DEBUG) {
            write_plate (image, ncols, nrows, p, out_path, out_image_name, i);
         }
      }
   }
   clock_t end = clock();
   printf("Time: %f\n", ((double)(end-start)) / ((double)(CLOCKS_PER_SEC / 1000.0)));
   return plates;
}
