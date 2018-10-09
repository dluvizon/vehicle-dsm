#ifndef _THOG_H_
#define _THOG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>

#include "utils.h"
#include "image.h"
#include "lanczos.h"
#include "clanczos.h"
#include "svm/svm.h"

#define TRUE 1
#define FALSE 0

#define L1 1

  typedef struct _struct_thog {
     int nh;                    /*Height image normalization (pixels).*/
     int ncx;                   /*Number of cells in x direction.*/
     int ncy;                   /*Number of cells in y direction.*/
     int noc;                   /*Total number of cells.*/
     int bpc;                   /*Bins per cell.*/
     int nob;                   /*Total number of bins.*/
     int norm;                  /*Image normalization.*/
     char wnorm[512];           /*Image normalization weight.*/
     double rad;                /*Image normalization weight radius*/
     char grad[512];            /*Image gradient option*/
     char hmetric[512];         /*Histogram normalization metric.*/
     char weight_function[512];
     int deformable_weights;
     int debug;
  } struct_thog;

  struct_thog load_settings (const char *settings);

double classify(const unsigned char *image, int nrows, int ncols,
		int x, int y, int w, int h, struct svm_model* model,
		double *prob, struct_thog sthog);

  double* thog (unsigned char *image, int nrows, int ncols, struct_thog sthog);

  void get_bin_pos (int bins_per_cell, double dtheta, int *bin,  double *factor);

  void get_mag_theta (double *grel, double *dnorm, double *dtheta, int rwidth, int rheight, int margin, double noise, char *gradient_option, int debug);

  void gradient_sobel (double *image, int width, int height, int x, int y, double *grad);
 
  void gradient_simple (double *image, int width, int height, int x, int y, double *grad);

  double cell_weight (char *weight_function, int ncz, int cz, int z, double zmax, double zmin);

  double StepFunc (int n, int k, double z);

  double Bernstein (int n, int k, double z);

  double BernsteinPoly (int n, int k, double z);

  double EdgeCore (int n, int k, double z);

  double Exp (int n, int k, double z);

  double gaussian (double z, double mu, double sigma);

  int choose_gaussian_weight_size (double dev);

  void compute_binomial_weights (double *weight, int rwt);

  void compute_gaussian_weights (double *weight, int rwt);

  void convert_to_log_scale (double *grey, int w, int h, double eps);

  double* get_gray_image (unsigned char *image, int w, int h);

  double *normalize_grey_image (double *grey, int w, int h, double *x_weight, int x_rwt, double *y_weight, int y_rwt, double noise);

  double get_grey_avg (double *grey, int w, int h, int x, int y, double *x_weight, int x_rwt, double *y_weight, int y_rwt);

  double get_grey_dev (double *grey, int w, int h, int x, int y, double *x_weight, int x_rwt, double *y_weight, int y_rwt, double AVG, double noise);

#endif
