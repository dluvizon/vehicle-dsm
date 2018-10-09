/**
 * @file detection-minetto.cpp
 * @author Rodrigo Minetto
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 10/03/2015
 */

#define _DLEVEL SNOOPERTEXT_DLEVEL

#include <assert.h>
#include <limits.h>
#include <math.h>

#include "detection-minetto.h"
#include "grouping/grouping.h"
#include "seg/toggle.h"
#include "config.h"
#include "lanes.h"
#include "debug.h"

#define VGRAY 255  /*GRAYLEVEL FOR VEHICLES*/
#define MGRAY 125  /*GRAYLEVEL FOR MOTORCYCLES*/
#define DEBUG 0

typedef struct _candidate {
   int xmin;
   int ymin;
   int xmax;
   int ymax;
   int npos;
   int total;
} candidate;

typedef struct _component {
   int xmin;
   int ymin;
   int ymax;
   int center;
   int group;
   bool motorcycle;
} component;

typedef struct _ncomponents {
   int size;
   int ngroups;
   component *array;
} ncomponents;

/*IO routines*/
static void save_frame_window_to_pgm_copy (
  uint8_t *img, 
  int imgw, 
  int w, int h, 
  const char *filename )
{
   int pix_cnt = 0;
   int mi = 0;
   FILE *file = fopen(filename, "w");
   if (!file) {
      print_err("could not open file %s", filename);
      return;
   }
   fprintf(file, "P2\n");
   fprintf(file, "%d %d\n", w, h);
   fprintf(file, "255\n");
   for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
         fprintf(file, "%d ", (int) img[mi]);
         if (pix_cnt % 12 == 0) {
            fprintf(file, "\n");
         }
         pix_cnt++;
         mi++;
      }
      mi += (imgw - w);
   }
   fclose(file);
}

static void write_components (
   int nrows, int ncols, 
   Vector cc,
   char *out_path, 
   char *out_image_name ) {

   unsigned char *out = (unsigned char *) malloc(nrows * ncols * sizeof(unsigned char));

   memset(out, 0, nrows * ncols * sizeof(unsigned char));

   for (int i = 0; i < vector_count(&cc); i++) {

      region *r = (region *) vector_get(&cc, i);

      for (int k = 0; k < r->np; k++) {
         int xx = (long) vector_get(&(r->xp), k);
         int yy = (long) vector_get(&(r->yp), k);
         out[yy * ncols + xx] = 125;
      }
      for (int k = r->box[0][0]; k <= r->box[0][1]; k++) {
         out[((int)r->center[1]) * ncols + k] = 255;
      }
      for (int k = r->box[0][0]; k <= r->box[0][1]; k++) {
         out[r->box[1][0] * ncols + k] = 255;
      }
      for (int k = r->box[0][0]; k <= r->box[0][1]; k++) {
         out[r->box[1][1] * ncols + k] = 255;
      }
   }
   write_uchar_to_pgm(out, ncols, nrows, "%s/%s_filtered.pgm", out_path, out_image_name);
   free(out);
}

static void write_chains (
   Vector** chains, 
   int nsets, 
   int wmin, int hmin,
   unsigned char *image, 
   int nrows, int ncols,
   char *out_path, 
   char *out_image_name )
{
   unsigned char *out = copy_uchar(image, nrows, ncols);
   for (int i = 0; i < nsets; i++) {
      int xmin = INT_MAX;
      int ymin = INT_MAX;
      int xmax = INT_MIN;
      int ymax = INT_MIN;
      /** Computing the chain dimensions based on the regions inside it: */
      for (int j = 0; j < vector_count(chains[i]); j++) {
         region *r = (region *)vector_get(chains[i], j);
         if (r->box[0][0] < xmin) { xmin = r->box[0][0]; }
         if (r->box[1][0] < ymin) { ymin = r->box[1][0]; }
         if (r->box[0][1] > xmax) { xmax = r->box[0][1]; }
         if (r->box[1][1] > ymax) { ymax = r->box[1][1]; }
      }
      int wbox = xmax - xmin;
      int hbox = ymax - ymin;
      if (wbox < wmin) { continue; }
      if (hbox < hmin) { continue; }
      if (hbox > wbox) { continue; }

      for (int j = 0; j < vector_count(chains[i]); j++) {
         region *r = (region *)vector_get(chains[i], j);
         for (int k = r->box[0][0]; k <= r->box[0][1]; k++) {
            out[r->box[1][0] * ncols + k] = 255;
            out[(r->box[1][0]-1) * ncols + k] = 255;
         }
         for (int k = r->box[0][0]; k <= r->box[0][1]; k++) {
            out[r->box[1][1] * ncols + k] = 255;
            out[(r->box[1][1]+1) * ncols + k] = 255;
         }
         for (int k = r->box[1][0]; k <= r->box[1][1]; k++) {
            out[k * ncols + r->box[0][0]] = 255;
            out[k * ncols + r->box[0][0]-1] = 255;
         }
         for (int k = r->box[1][0]; k <= r->box[1][1]; k++) {
            out[k * ncols + r->box[0][1]] = 255;
            out[k * ncols + r->box[0][1]+1] = 255;
         }
      }
   }
   write_uchar_to_pgm(out, ncols, nrows, "%s/%s_chains.pgm", out_path, out_image_name);
   free(out);
}

/*End IO routines!*/

void free_vector_elemets (Vector *v) {
   while (v->count--) {
      region *r = (region *) v->data[v->count];
      vector_free(&(r->xp));
      vector_free(&(r->yp));
      free(r);
   }
}

void Sobel (int xmin, int ymin, int xmax, int ymax, unsigned char *image, int nrows, int ncols, unsigned char *out, int strenght, int GRAY) {

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

   for (int y = ymin; y < ymax; y++) {
      for (int x = xmin; x < xmax; x++) {
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

   for (int y = ymin; y < ymax; y++) {
      for (int x = xmin; x < xmax; x++) {
         if (tmp[y * ncols + x] < threshold) {
            out[y * ncols + x] = 0;
         }
         else {
            out[y * ncols + x] = GRAY;
         }
      }
   }

   free(tmp);
}

void Edge_Filtering (int xmin, int ymin, int xmax, int ymax, int nrows, int ncols, unsigned char *out, double rmin, double rmax, int GRAY) {

   /*Excluding every component that touches in the edges!*/
   for (int x = xmin; x <= xmax; x++) {
      out[ymin * ncols + x] = 255;
      out[ymax * ncols + x] = 255;
   }
   for (int y = ymin; y <= ymax; y++) {
      out[y * ncols + xmin] = 255;
      out[y * ncols + xmax] = 255;
   }
   /*End*/

   Vector cc = find_connected_components4 (out, nrows, ncols, rmin, rmax);

   for (int y = ymin; y <= ymax; y++) {
      for (int x = xmin; x <= xmax; x++) {
 	out[y * ncols + x] = 0;
      }
   }
   for (int i = 0; i < vector_count(&cc); i++) {
      region *r = (region *) vector_get(&cc, i);
      for (int k = 0; k < r->np; k++) {
         int xx = (long) vector_get(&(r->xp), k);
	 int yy = (long) vector_get(&(r->yp), k);
	 out[yy * ncols + xx] = GRAY;
      }
   }
   free_vector_elemets(&cc);
   vector_free(&cc);
}


/**/
static void Dilate (int xmin, int ymin, int xmax, int ymax, unsigned char *image, int nrows, int ncols, unsigned char *out, int mask[], int mwidth, int GRAY) {

   xmin = max(xmin, mwidth/2);
   ymin = max(ymin, 0);
   xmax = min(xmax, ncols - mwidth/2);
   ymax = min(ymax, nrows);

   for (int y = ymin; y < ymax; y++) {
      for (int x = xmin; x < xmax; x++) {
         int pimg = y * ncols + x;
         if (image[pimg] > 0) {
            for (int i = -mwidth/2; i <= +mwidth/2; i++) {
               int pout = y * ncols + (x + i);
               out[pout] = GRAY;
            }
         }
      }
   }
}

/**/
static ncomponents* Find_Regions (
   Vector** chains, 
   int nsets,
   int wmin, int hmin, 
   int step, 
   int nrows, int ncols, 
   unsigned char *image ) {

   bool motorcyle = false;
   ncomponents *c = (ncomponents *) malloc(sizeof(ncomponents));
   c->size = 0;
   c->ngroups = 0;
   c->array = (component*) malloc((nrows * ncols) * sizeof(component));
   int group = 0;

   unsigned char *untouch = (unsigned char *)malloc((nrows * ncols) * sizeof(unsigned char));
   memset (untouch, 0, nrows * ncols * sizeof(unsigned char));

   for (int i = 0; i < nsets; i++) {
      int xmin = INT_MAX;
      int ymin = INT_MAX;
      int xmax = INT_MIN;
      int ymax = INT_MIN;
      /*Computing the chain dimensions based on the regions inside it: */
      for (int j = 0; j < vector_count(chains[i]); j++) {
         region *r = (region *) vector_get(chains[i], j);
         if (r->box[0][0] < xmin) { xmin = r->box[0][0]; }
         if (r->box[1][0] < ymin) { ymin = r->box[1][0]; }
         if (r->box[0][1] > xmax) { xmax = r->box[0][1]; }
         if (r->box[1][1] > ymax) { ymax = r->box[1][1]; }
         for (int k = 0; k < r->np; k++) {
            int xx = (long) vector_get(&(r->xp), k);
            int yy = (long) vector_get(&(r->yp), k);
            untouch[yy * ncols + xx] = 255;
            if (image[yy * ncols + xx] == MGRAY) {
               motorcyle = true;
            }
         }
      }
      int wbox = xmax - xmin;
      int hbox = ymax - ymin;
      if (wbox < wmin) { continue; }
      if (hbox < hmin) { continue; }
      if (hbox > wbox) { continue; }

      for (int j = 0; j < vector_count(chains[i]); j++) {
         region *r = (region *) vector_get(chains[i], j);
         if (motorcyle) {
            for (int k = 0; k < r->np; k++) {
               int xx = (long) vector_get(&(r->xp), k);
               int yy = (long) vector_get(&(r->yp), k);
               c->array[c->size].xmin = xx; 
               c->array[c->size].center = yy;
               c->array[c->size].group = group;
               c->array[c->size].motorcycle = true;
               c->size++;
            }
         }
         else {
           for (int k = r->box[0][0]; k <= r->box[0][1]; k += step) {
              c->array[c->size].xmin = k;
              c->array[c->size].ymin = r->box[1][0];
              c->array[c->size].ymax = r->box[1][1];
               int min = c->array[c->size].ymax, max = c->array[c->size].ymin;
               for (int t = c->array[c->size].ymin; t < c->array[c->size].ymax; t++) {
                  int pixel = t * ncols + k;
                  if (untouch[pixel] > 0) {
                     if (t < min) { min = t; }
                     if (t > max) { max = t; }
                  }
               }
               c->array[c->size].center = min + (max - min)/2;
               c->array[c->size].group = group;
               c->array[c->size].motorcycle = false;
               c->size++;
            }
         }
      }
      group++;
    }
    c->ngroups = group;
    free (untouch);
    return c;
}

static vector<plate> THOG_Classification (
   unsigned char *image,
   int nrows, int ncols, 
   ncomponents *c, 
   struct svm_model* model,
   double *prob, 
   struct_thog *sthog, 
   char *out_path, 
   char *out_image_name, 
   int npoints, 
   int wbox, int hbox ) {

   candidate *samples = (candidate *)malloc(c->ngroups * sizeof(candidate));

   for (int i = 0; i < c->ngroups; i++) {
      samples[i].xmin = 99999;
      samples[i].xmax = 0;
      samples[i].ymin = 99999;
      samples[i].ymax = 0;
      samples[i].npos = 0;
      samples[i].total = 0;
   }

   unsigned char *map = NULL;
   if (DEBUG) {
     map = (unsigned char *)malloc(nrows * ncols * sizeof(unsigned char));
     memset (map, 0, nrows * ncols * sizeof(unsigned char));
   }

   int n = 0;
   while (n < c->size) {
     if (c->array[n].motorcycle) { wbox = 24; hbox = 12; }
     else { wbox = 48; hbox = 24; }
     if ( (c->array[n].xmin - wbox/2) > 0 && 
          (c->array[n].center - hbox/2) > 0 && 
          (c->array[n].center + hbox/2) < nrows && 
          (c->array[n].xmin + wbox/2) < ncols ) {

         if ( (classify (image, nrows, ncols, c->array[n].xmin - wbox/2, c->array[n].center - hbox/2, wbox, hbox, model, prob, *sthog) > 0.0) ) {
				
            samples[c->array[n].group].npos++;
            if (samples[c->array[n].group].xmin > (c->array[n].xmin   - wbox/2)) {
               samples[c->array[n].group].xmin = c->array[n].xmin   - wbox/2;
	    }
	    if (samples[c->array[n].group].ymin > (c->array[n].center - hbox/2)) {
               samples[c->array[n].group].ymin = c->array[n].center - hbox/2;
            }
            if (samples[c->array[n].group].xmax < (c->array[n].xmin   + wbox/2)) {
               samples[c->array[n].group].xmax = c->array[n].xmin   + wbox/2;
            }
            if (samples[c->array[n].group].ymax < (c->array[n].center + hbox/2)) {
               samples[c->array[n].group].ymax = c->array[n].center + hbox/2;
            }
            if (DEBUG) {
               map[ncols*c->array[n].ymin + c->array[n].xmin] = 255;
               map[ncols*c->array[n].center + c->array[n].xmin] = 255;
            }
         } else {
            if (DEBUG) {
               map[ncols*c->array[n].ymin + c->array[n].xmin] = 100;
               map[ncols*c->array[n].center + c->array[n].xmin] = 100;
            }
         }
      }
      samples[c->array[n].group].total++;
      n++;
   }

   vector<plate> plates;
   for (int i = 0; i < c->ngroups; i++) {
      if (samples[i].npos > npoints) {
         plate p = {samples[i].xmin, samples[i].ymin, (samples[i].xmax - samples[i].xmin + 1), (samples[i].ymax - samples[i].ymin + 1) };
         p.speed = 0.0;
         p.from = NULL;
         plates.push_back(p);
      }
   }
      
   if (DEBUG) {
      unsigned char *copy = copy_uchar(image, nrows, ncols);
      char out_filename[512];
      sprintf(out_filename, "%s/%s_detection.txt", out_path, out_image_name);
      FILE *file = fopen(out_filename, "w");
      for (int i = 0; i < c->ngroups; i++) {
         if (samples[i].npos > npoints) {
            for (int k = samples[i].xmin; k < samples[i].xmax; k++) {
               copy[ncols*samples[i].ymin + k] = 255;
               copy[ncols*samples[i].ymax + k] = 255;
               copy[ncols*(samples[i].ymin-1) + k] = 0;
               copy[ncols*(samples[i].ymax+1) + k] = 0;
            }
            for (int k = samples[i].ymin; k < samples[i].ymax; k++) {
               copy[ncols*k + samples[i].xmin] = 255;
               copy[ncols*k + samples[i].xmax] = 255;
               copy[ncols*k + samples[i].xmin - 1] = 0;
               copy[ncols*k + samples[i].xmax + 1] = 0;
            }
         }
      }
      fclose (file);
      write_uchar_to_pgm (copy, ncols, nrows, "%s/%s_detection.pgm", out_path, out_image_name);
      write_uchar_to_pgm (map, ncols, nrows, "%s/%s_map.pgm", out_path, out_image_name);
      free(map);
      free(copy);
   }
   free(samples);
   return plates;
}

vector<plate> detection_minetto(
   unsigned char *image,
   struct snoopertext *s,
   int ncols, int nrows,
   char *out_path,
   char *out_image_name,
   vector<struct sub_slope> &slopes,
   int margin_X, int margin_Y,
   int iframe,
   unsigned char *umask) {

   clock_t start = clock();

   /*Parameter settings:*/
      /*Sobel*/
      int mstrenght = 4;                     /*change*/
      int vstrenght = 2;                     /*change*/
      int mwidth = 7;                        /*change*/  
      int mask[7] = {1, 1, 1, 1, 1, 1, 1};   /*change*/
      /*Character filtering*/
      double mchar_rmin = 4;
      double mchar_rmax = 300;
      double vchar_rmin = 4;
      double vchar_rmax = 120;
      /*Grouping*/
      double t1 = 0.7;                       /*fixed*/
      double t2 = 1.1;                       /*fixed*/
      double t3 = 0.4;                       /*fixed*/
      /*Grouping filtering*/
      int wmin = 32;                         /*fixed*/
      int hmin = 10;                         /*fixed*/
      double reg_rmin = 10;                  /*fixed*/
      double reg_rmax = 200;                 /*fixed*/ 
      /*T-HOG*/
      int window_w = 48;                     /*fixed*/
      int window_h = 24;                     /*fixed*/
      int step = 4;                          /*fixed*/ 
      int npoints = 8;                       /*change*/
   /*End settings*/

   vector<plate> plates;

   struct svm_model* model = s->model;

   double *prob = s->prob;

   struct_thog *sthog = &s->sthog;

   int cter = 0;
   for (int i = 0; i < slopes.size(); i++) {
      if (!slopes.at(i).search_plate) { continue; }
      cter++;
   }
   if (cter == 0) { return plates; }

   memset(s->edges, 0, nrows * ncols * sizeof(unsigned char));

   memset(s->dilate, 0, nrows * ncols * sizeof(unsigned char));
  
   for (int i = 0; i < slopes.size(); i++) {

      if (!slopes.at(i).search_plate) { continue; }

      if (DEBUG) {
         char filename[256];
         sprintf(filename, "%s/%05d_s%02d.pgm", out_path, iframe, i);
         save_frame_window_to_pgm_copy (image + slopes.at(i).left + (slopes.at(i).up * ncols), ncols, slopes.at(i).right - slopes.at(i).left, slopes.at(i).down - slopes.at(i).up, filename);
      }

      if ( (slope_is_moto(slopes.at(i))) && (slopes.at(i).ntries > 2) ) {
         Sobel (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, image, nrows, ncols, s->edges, mstrenght, MGRAY);
      } 
      else {
         Sobel (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, image, nrows, ncols, s->edges, vstrenght, VGRAY);
      } 

      if (DEBUG) {
         write_uchar_to_pgm (s->edges, ncols, nrows, "%s/%s_%02d_sobel.pgm", out_path, out_image_name, i);
      }
                
      if ( (slope_is_moto(slopes.at(i))) && (slopes.at(i).ntries > 2) ) {
         Edge_Filtering (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, nrows, ncols, s->edges, mchar_rmin, mchar_rmax, MGRAY);
      }
      else {
         Edge_Filtering (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, nrows, ncols, s->edges, vchar_rmin, vchar_rmax, VGRAY);
      }

      if (DEBUG) {
         write_uchar_to_pgm (s->edges, ncols, nrows, "%s/%s_%02d_filtered.pgm", out_path, out_image_name, i);
      }

      if ( (slope_is_moto(slopes.at(i))) && (slopes.at(i).ntries > 2) ) {
         Dilate (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, s->edges, nrows, ncols, s->dilate, mask, mwidth, MGRAY);
      }
      else {
         Dilate (slopes.at(i).left, slopes.at(i).up, slopes.at(i).right, slopes.at(i).down, s->edges, nrows, ncols, s->dilate, mask, mwidth, VGRAY);
      }

      if (DEBUG) {
         write_uchar_to_pgm(s->dilate, ncols, nrows, "%s/%s_%02d_dilate.pgm", out_path, out_image_name, i);
      }
   }

   Vector ncc = find_connected_components4 (s->dilate, nrows, ncols, reg_rmin, reg_rmax);

   int nsets;
   Vector** nchains = make_chains (ncc, t1, t2, t3, &nsets);

   ncomponents *comp = Find_Regions (nchains, nsets, wmin, hmin, step, nrows, ncols, s->dilate);

   plates = THOG_Classification (image, nrows, ncols, comp, model, prob, sthog, out_path, out_image_name, npoints, window_w, window_h);

   if (DEBUG) {
      write_components (nrows, ncols, ncc, out_path, out_image_name);
      write_chains (nchains, nsets, wmin, hmin, image, nrows, ncols, out_path, out_image_name);
   }

   free_vector_elemets(&ncc);
   vector_free(&ncc);
   for (int i = 0; i < nsets; i++) {
      vector_free(nchains[i]);
   }
   free(nchains);
   free(comp->array);
   free(comp);

   clock_t end = clock();

   printf("Time: %f\n", ((double)(end-start))/((double)(CLOCKS_PER_SEC/1000.0)));

   return plates;
}

