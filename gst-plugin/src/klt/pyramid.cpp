/**
 * pyramid.c
 *
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include "base.h"
#include "error.h"
#include "convolve.h"
#include "pyramid.h"

static void print_double(double num)
{
	if (num < 10.) {
		printf("   ");
	} else if (num < 100.) {
		printf("  ");
	} else if (num < 1000.) {
		printf(" ");
	}
	printf("%.3f", num);
}

static clock_t _clk0;
static clock_t _clk1;

static void chrono_reset()
{
	_clk0 = clock();
}

static void chrono_show(char *msg)
{
	double msdt;

	_clk1 = clock();
	msdt = (_clk1 - _clk0) * 1000. / CLOCKS_PER_SEC;
	print_double(msdt);
	printf(" -- %s\n", msg);
}

static void chrono_show_and_reset(char *msg)
{
	double msdt;

	_clk1 = clock();
	msdt = (_clk1 - _clk0) * 1000. / CLOCKS_PER_SEC;
	chrono_reset();
	print_double(msdt);
	printf(" -- %s\n", msg);
}

/*********************************************************************
 *
 */

_KLT_Pyramid _KLTCreatePyramid(
  int ncols,
  int nrows,
  int subsampling,
  int nlevels)
{
  _KLT_Pyramid pyramid;
  int nbytes = sizeof(_KLT_PyramidRec) +	
    nlevels * sizeof(_KLT_FloatImage *) +
    nlevels * sizeof(int) +
    nlevels * sizeof(int);
  int i;

  if (subsampling != 2 && subsampling != 4 && 
      subsampling != 8 && subsampling != 16 && subsampling != 32)
    KLTError("(_KLTCreatePyramid)  Pyramid's subsampling must "
             "be either 2, 4, 8, 16, or 32");

     
  /* Allocate memory for structure and set parameters */
 // printf("%s: MALLOC\n", __func__);
  pyramid = (_KLT_Pyramid)  malloc(nbytes);
  if (pyramid == NULL) {
    KLTError("(_KLTCreatePyramid)  Out of memory");
  }
     
  /* Set parameters */
  pyramid->subsampling = subsampling;
  pyramid->nLevels = nlevels;
  pyramid->img = (_KLT_FloatImage *) (pyramid + 1);
  pyramid->ncols = (int *) (pyramid->img + nlevels);
  pyramid->nrows = (int *) (pyramid->ncols + nlevels);

  /* Allocate memory for each level of pyramid and assign pointers */
  for (i = 0 ; i < nlevels ; i++)  {
    pyramid->img[i] =  _KLTCreateFloatImage(ncols, nrows);
    pyramid->ncols[i] = ncols;  pyramid->nrows[i] = nrows;
    ncols /= subsampling;  nrows /= subsampling;
  }

  return pyramid;
}


/*********************************************************************
 *
 */

void _KLTFreePyramid(
  _KLT_Pyramid pyramid)
{
  int i;

  /* Free images */
  for (i = 0 ; i < pyramid->nLevels ; i++)
    _KLTFreeFloatImage(pyramid->img[i]);

  /* Free structure */
  free(pyramid);
}

/*********************************************************************
 *
 */

_KLT_Pyramid _KLTCopyPyramid (
  _KLT_Pyramid src )
{
  int i;

  _KLT_Pyramid dst = _KLTCreatePyramid (src->ncols[0], src->nrows[0], src->subsampling, src->nLevels);

  /* Free images */
  for (i = 0 ; i < src->nLevels ; i++) {
    _KLTCopyFloatImage (src->img[i], dst->img[i]);
  }  

  return dst;
}

/*********************************************************************
 *
 */

void _KLTComputePyramid(_KLT_FloatImage img, _KLT_Pyramid pyramid,
		float sigma_fact)
{
	_KLT_FloatImage currimg;
	_KLT_FloatImage tmpimg;
	int ncols = img->ncols;
	int nrows = img->nrows;
	int subsampling = pyramid->subsampling;
	int subhalf = subsampling / 2;
	float sigma = subsampling * sigma_fact;	/* empirically determined */
	int oldncols;
	int i, x, y;

	if (subsampling != 2 && subsampling != 4 && subsampling != 8 &&
			subsampling != 16 && subsampling != 32) {
		KLTError("(_KLTComputePyramid)  Pyramid's subsampling must "
				"be either 2, 4, 8, 16, or 32");
	}

	assert(pyramid->ncols[0] == img->ncols);
	assert(pyramid->nrows[0] == img->nrows);

	/* Copy original image to level 0 of pyramid */
	memcpy(pyramid->img[0]->data, img->data,
			ncols * nrows * sizeof(float));

	currimg = img;
	for (i = 1 ; i < pyramid->nLevels ; i++) {

		tmpimg = _KLTCreateFloatImage(ncols, nrows);

		//_KLTWriteFloatImageToPGM(currimg, "saida/before_smooth.pgm");

		_KLTComputeSmoothedImage(currimg, sigma, tmpimg);

		//_KLTWriteFloatImageToPGM(tmpimg, "saida/after_smooth.pgm");

		/* Subsample */
		oldncols = ncols;
		ncols /= subsampling;
		nrows /= subsampling;
		for (y = 0 ; y < nrows ; y++) {
			for (x = 0 ; x < ncols ; x++) {
				pyramid->img[i]->data[y * ncols+x] =
					tmpimg->data[(subsampling * y +
							subhalf) * oldncols +
					(subsampling * x + subhalf)];
			}
		}

		/* Reassign current image */
		currimg = pyramid->img[i];
		_KLTFreeFloatImage(tmpimg);
	}
}
 
/*********************************************************************
 *
 */

void _KLTComputePyramidC (
  vector<KLT_FeatureList> featurelist, 
  _KLT_FloatImage img, 
  _KLT_Pyramid pyramid,
  float sigma_fact,
  KLT_TrackingContext tc)
{
  _KLT_FloatImage currimg, tmpimg;
  int ncols = img->ncols, nrows = img->nrows;
  int subsampling = pyramid->subsampling;
  int subhalf = subsampling / 2;
  float sigma = subsampling * sigma_fact;  /* empirically determined */
  int oldncols;
  int i, x, y;
	
  if (subsampling != 2 && subsampling != 4 && 
      subsampling != 8 && subsampling != 16 && subsampling != 32)
    KLTError("(_KLTComputePyramid)  Pyramid's subsampling must "
             "be either 2, 4, 8, 16, or 32");

  assert(pyramid->ncols[0] == img->ncols);
  assert(pyramid->nrows[0] == img->nrows);

  /* Copy original image to level 0 of pyramid */
  memcpy(pyramid->img[0]->data, img->data, ncols*nrows*sizeof(float));

  currimg = img;
  for (i = 1 ; i < pyramid->nLevels ; i++)  {
    tmpimg = _KLTCreateFloatImage(ncols, nrows);
    _KLTComputeSmoothedImageC (featurelist, currimg, sigma, tmpimg, i, tc);


    /* Subsample */
    oldncols = ncols;
    ncols /= subsampling;  nrows /= subsampling;
    for (y = 0 ; y < nrows ; y++)
      for (x = 0 ; x < ncols ; x++)
        pyramid->img[i]->data[y*ncols+x] = 
          tmpimg->data[(subsampling*y+subhalf)*oldncols +
                      (subsampling*x+subhalf)];

    /* Reassign current image */
    currimg = pyramid->img[i];
				
    _KLTFreeFloatImage(tmpimg);
  }
}











