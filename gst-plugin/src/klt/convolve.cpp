/**
 * convolve.c
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "base.h"
#include "error.h"
#include "convolve.h"
#include "klt_util.h"

#define MAX_KERNEL_WIDTH 71

typedef struct {
	int width;
	float data[MAX_KERNEL_WIDTH];
} ConvolutionKernel;

/* Kernels */
static ConvolutionKernel gauss_kernel;
static ConvolutionKernel gaussderiv_kernel;
static float sigma_last = -10.0;

/**
 * _KLTToFloatImage
 * Given a pointer to image data (probably unsigned chars), copy
 * data to a float image.
 */
void _KLTToFloatImage(const KLT_PixelType *img, int ncols, int nrows,
		_KLT_FloatImage floatimg)
{
	const KLT_PixelType *ptrend = img + ncols*nrows;
	float *ptrout = floatimg->data;

	/* Output image must be large enough to hold result. */
	assert(floatimg->ncols >= ncols);
	assert(floatimg->nrows >= nrows);

	floatimg->ncols = ncols;
	floatimg->nrows = nrows;

	while (img < ptrend) {
		*ptrout++ = (float) *img++;
	}
}


/**
 * _computeKernels
 */
static void _computeKernels(float sigma, ConvolutionKernel *gauss,
		ConvolutionKernel *gaussderiv)
{
	const float factor = 0.01f;	/* For truncating tail. */
	int i;

	assert(MAX_KERNEL_WIDTH % 2 == 1);
	assert(sigma >= 0.0);

	/* Compute kernels, and automatically determine widths */
	{
		const int hw = MAX_KERNEL_WIDTH / 2;
		float max_gauss = 1.0f, max_gaussderiv =
			(float) (sigma*exp(-0.5f));

		/* Compute gauss and deriv. */
		for (i = -hw ; i <= hw ; i++) {
			gauss->data[i + hw] =
				(float) exp(-i * i / (2 * sigma * sigma));
			gaussderiv->data[i + hw] = -i * gauss->data[i + hw];
		}

		/* Compute widths. */
		gauss->width = MAX_KERNEL_WIDTH;
		for (i = -hw; fabs(gauss->data[i + hw] / max_gauss) < factor;
				i++, gauss->width -= 2);
		gaussderiv->width = MAX_KERNEL_WIDTH;

		for (i = -hw; fabs(gaussderiv->data[i + hw] / max_gaussderiv)
				< factor; i++, gaussderiv->width -= 2);

		if (gauss->width == MAX_KERNEL_WIDTH ||
				gaussderiv->width == MAX_KERNEL_WIDTH) {
			KLTError("(_computeKernels) MAX_KERNEL_WIDTH %d is "
					"too small for a sigma of %f",
					MAX_KERNEL_WIDTH, sigma);
		}
	}

	/* Shift if width less than MAX_KERNEL_WIDTH. */
	for (i = 0; i < gauss->width; i++) {
		gauss->data[i] = gauss->data[i +
			(MAX_KERNEL_WIDTH - gauss->width) / 2];
	}

	for (i = 0; i < gaussderiv->width; i++) {
		gaussderiv->data[i] = gaussderiv->data[i +
			(MAX_KERNEL_WIDTH - gaussderiv->width) / 2];
	}
	/* Normalize gauss and deriv. */
	{
		const int hw = gaussderiv->width / 2;
		float den;

		den = 0.0;
		for (i = 0 ; i < gauss->width; i++) {
			den += gauss->data[i];
		}
		for (i = 0; i < gauss->width; i++) {
			gauss->data[i] /= den;
		}
		den = 0.0;
		for (i = -hw; i <= hw; i++) {
			den -= i*gaussderiv->data[i + hw];
		}
		for (i = -hw; i <= hw; i++) {
			gaussderiv->data[i + hw] /= den;
		}
	}

	sigma_last = sigma;
}
	
/**
 * _KLTGetKernelWidths
 */
void _KLTGetKernelWidths(float sigma, int *gauss_width, int *gaussderiv_width)
{
	_computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
	*gauss_width = gauss_kernel.width;
	*gaussderiv_width = gaussderiv_kernel.width;
}

/**
 * _convolveImageHoriz
 */
static void _convolveImageHoriz(_KLT_FloatImage imgin,
		ConvolutionKernel kernel, _KLT_FloatImage imgout)
{
	float *ptrrow = imgin->data; /* Points to row's first pixel */
	float *ptrout; /* Points to next output pixel */
	float *ppp;
	float sum;
	int radius = kernel.width / 2;
	int ncols = imgin->ncols;
	int nrows = imgin->nrows;
	register int i, j, k;

	ptrout = imgout->data;

	/* Kernel width must be odd */
	assert(kernel.width % 2 == 1);

	/* Must read from and write to different images */
	assert(imgin != imgout);

	/* Output image must be large enough to hold result */
	assert(imgout->ncols >= imgin->ncols);
	assert(imgout->nrows >= imgin->nrows);

	/* For each row, do ... */
	for (j = 0 ; j < nrows ; j++)  {
		/* Zero leftmost columns */
		for (i = 0 ; i < radius ; i++) {
			*ptrout++ = 0.0;
		}

		/* Convolve middle columns with kernel */
		for (; i < ncols - radius ; i++)  {
			ppp = ptrrow + i - radius;
			sum = 0.0;
			for (k = kernel.width - 1; k >= 0; k--) {
				sum += *ppp++ * kernel.data[k];
			}
			*ptrout++ = sum;
		}

		/* Zero rightmost columns */
		for ( ; i < ncols ; i++) {
			*ptrout++ = 0.0;
		}
		ptrrow += ncols;
	}
}

/**
 * _convolveImageHoriz
 */
static void _convolveImageHorizB (
  _KLT_FloatImage imgin,
  ConvolutionKernel kernel,
  _KLT_FloatImage imgout,
  plate pl)
{
  float *ptrrow = imgin->data;           /* Points to row's first pixel */
  register float *ptrout = imgout->data, /* Points to next output pixel */
    *ppp;
  register float sum;
  register int radius = kernel.width / 2;
  register int ncols = imgin->ncols, nrows = imgin->nrows;
  register int i, j, k;

  /* Kernel width must be odd */
  assert(kernel.width % 2 == 1);

  /* Must read from and write to different images */
  assert(imgin != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= imgin->ncols);
  assert(imgout->nrows >= imgin->nrows);

  unsigned int p;

  /*TO-DO: CUIDADO colocar IF para não ultrapassar devido ao radius */

  /* For each plate, do ... */
  //for (p = 0; p < plates.size(); p++) {
    /* For each plate row, do ... */
    for (j = pl.y; j < pl.y + pl.h; j++) {
      /* Convolve columns with kernel */
      for (i = pl.x; i < pl.x + pl.w; i++)  {

        if (
          (j * ncols + i - radius < 0) ||
          (j * ncols + i - radius + kernel.width) >= (ncols * nrows)
        ) { continue; }

        ppp = ptrrow + j * ncols + i - radius;
        sum = 0.0;
        for (k = kernel.width-1 ; k >= 0 ; k--)
          sum += *ppp++ * kernel.data[k];
        *(ptrout + j * ncols + i) = sum;
      }
    } 
  //printf("Depois da iteração _convolveImageHorizB\n");
  //}
}

/*********************************************************************
 * _convolveImageHoriz
 */
static void _convolveImageHorizC(vector<KLT_FeatureList> featurelist,
		_KLT_FloatImage imgin, ConvolutionKernel kernel,
		_KLT_FloatImage imgout, int level, KLT_TrackingContext tc)
{
	float *ptrrow = imgin->data;	/* Points to row's first pixel. */
	 /* Points to next output pixel. */
	register float *ptrout = imgout->data;
	register float *ppp;
	register float sum;
	register int radius = kernel.width / 2;
	register int ncols = imgin->ncols, nrows = imgin->nrows;
	register int i, j, k;

	/* Kernel width must be odd. */
	assert(kernel.width % 2 == 1);

	/* Must read from and write to different images. */
	assert(imgin != imgout);

	/* Output image must be large enough to hold result. */
	assert(imgout->ncols >= imgin->ncols);
	assert(imgout->nrows >= imgin->nrows);

	int s = imgin->ncols * imgin->nrows;
	unsigned int p;
	int margin = 50;

	/* For each plate, do ... */
	for (p = 0; p < featurelist.size(); p++) {
		int xmin = featurelist.at(p)->xmin / pow(tc->subsampling,
				level) + featurelist.at(p)->dx - margin;
		int xmax = featurelist.at(p)->xmax / pow(tc->subsampling,
				level) + featurelist.at(p)->dx + margin;
		int ymin = featurelist.at(p)->ymin / pow(tc->subsampling,
				level) + featurelist.at(p)->dy - margin;
		int ymax = featurelist.at(p)->ymax / pow(tc->subsampling,
				level) + featurelist.at(p)->dy + margin; 

		/* For each plate row, do ... */
		for (j = ymin; j < ymax; j++) {
			/* Convolve columns with kernel. */
			for (i = xmin; i < xmax; i++)  {
				int p1 = j * ncols + i - radius;
				int p2 = j * ncols + i + radius;
				if ((p1 >= 0) && (p2 < s)) {
					ppp = ptrrow + j * ncols + i - radius;
					sum = 0.0;
					for (k = kernel.width - 1; k >= 0;
							k--) {
						sum += *ppp++ * kernel.data[k];
					}
					*(ptrout + j * ncols + i) = sum;
				}
			}
		} 
	}
}

/*********************************************************************
 * _convolveImageVert
 */
static void _convolveImageVert(_KLT_FloatImage imgin,
		ConvolutionKernel kernel, _KLT_FloatImage imgout)
{
	float *ptrcol = imgin->data; /* Points to row's first pixel */
	float *ptrout = imgout->data; /* Points to next output pixel */
	float *ppp;
	float sum;
	int radius = kernel.width / 2;
	int ncols = imgin->ncols, nrows = imgin->nrows;
	register int i, j, k;

	/* Kernel width must be odd */
	assert(kernel.width % 2 == 1);

	/* Must read from and write to different images */
	assert(imgin != imgout);

	/* Output image must be large enough to hold result */
	assert(imgout->ncols >= imgin->ncols);
	assert(imgout->nrows >= imgin->nrows);

	/* For each column, do ... */
	for (i = 0 ; i < ncols ; i++) {
		/* Zero topmost rows */
		for (j = 0; j < radius; j++) {
			*ptrout = 0.0;
			ptrout += ncols;
		}

		/* Convolve middle rows with kernel */
		for (; j < nrows - radius; j++) {
			ppp = ptrcol + ncols * (j - radius);
			sum = 0.0;
			for (k = kernel.width - 1 ; k >= 0 ; k--) {
				sum += *ppp * kernel.data[k];
				ppp += ncols;
			}
			*ptrout = sum;
			ptrout += ncols;
		}

		/* Zero bottommost rows */
		for (; j < nrows; j++) {
			*ptrout = 0.0;
			ptrout += ncols;
		}

		ptrcol++;
		ptrout -= nrows * ncols - 1;
	}
}

/*********************************************************************
 * _convolveImageVert
 */

static void _convolveImageVertB (
  _KLT_FloatImage imgin,
  ConvolutionKernel kernel,
  _KLT_FloatImage imgout,
  plate pl)
{
  float *ptrcol = imgin->data;            /* Points to row's first pixel */
  register float *ptrout = imgout->data,  /* Points to next output pixel */
    *ppp;
  register float sum;
  register int radius = kernel.width / 2;
  register int ncols = imgin->ncols, nrows = imgin->nrows;
  register int i, j, k;

  /* Kernel width must be odd */
  assert(kernel.width % 2 == 1);

  /* Must read from and write to different images */
  assert(imgin != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= imgin->ncols);
  assert(imgout->nrows >= imgin->nrows);

  unsigned int p;

  /*TO-DO: CUIDADO colocar IF para não ultrapassar devido ao radius */

  /* For each plate, do ... */
  //for (p = 0; p < plates.size(); p++) {
    /* For each plate column, do ... */
    for (i = pl.x; i < pl.x + pl.w; i++)  {
      /* Convolve rows with kernel */
      for (j = pl.y; j < pl.y + pl.h; j++) {
        if (
          ((j - radius) * ncols + i < 0) ||
           (j - radius + kernel.width) * ncols + i >= (ncols * nrows)
        ) { continue; }
        //printf("depois: %d\n", (j - radius) * ncols + i);
        ppp = ptrcol + (j - radius) * ncols + i;
        //if ((j - radius) * ncols + i < 0) { printf("caiu fora\n"); }
        sum = 0.0;
        for (k = kernel.width-1 ; k >= 0 ; k--) {
          sum += *ppp * kernel.data[k];
          ppp += ncols;
        }
        *(ptrout + j * ncols + i) = sum;
      }
    } 
  //}
  /*TO-DO: CUIDADO colocar IF para não ultrapassar devido ao radius */
}

static void _convolveImageVertC(vector<KLT_FeatureList> featurelist,
		_KLT_FloatImage imgin, ConvolutionKernel kernel,
		_KLT_FloatImage imgout, int level, KLT_TrackingContext tc)
{
	float *ptrcol = imgin->data;	/* Points to row's first pixel. */
	float *ptrout = imgout->data;	/* Points to next output pixel. */
	float *ppp;
	register float sum;
	int radius = kernel.width / 2;
	int ncols = imgin->ncols;
	int nrows = imgin->nrows;
	register int i, j, k;

	/* Kernel width must be odd. */
	assert(kernel.width % 2 == 1);

	/* Must read from and write to different images. */
	assert(imgin != imgout);

	/* Output image must be large enough to hold result. */
	assert(imgout->ncols >= imgin->ncols);
	assert(imgout->nrows >= imgin->nrows);

	int s = imgin->ncols * imgin->nrows;
	int margin = 50;

	/* For each plate, do ... */
	for (unsigned int p = 0; p < featurelist.size(); p++) {
		int xmin = featurelist.at(p)->xmin / pow(tc->subsampling,
				level) + featurelist.at(p)->dx - margin;
		int xmax = featurelist.at(p)->xmax / pow(tc->subsampling,
				level) + featurelist.at(p)->dx + margin;
		int ymin = featurelist.at(p)->ymin / pow(tc->subsampling,
				level) + featurelist.at(p)->dy - margin;
		int ymax = featurelist.at(p)->ymax / pow(tc->subsampling,
				level) + featurelist.at(p)->dy + margin;

		/* For each plate column, do ... */
		for (i = xmin; i < xmax; i++)  {
			/* Convolve rows with kernel. */
			for (j = ymin; j < ymax; j++) {
				int p1 = (j - radius) * ncols + i;
				int p2 = (j + radius) * ncols + i;
				if ((p1 >= 0) && (p2 < s)) {
					ppp = ptrcol + (j - radius) * ncols + i;
					sum = 0.0;
					for (k = kernel.width - 1; k >= 0;
							k--) {
						sum += *ppp * kernel.data[k];
						ppp += ncols;
					}
					*(ptrout + j * ncols + i) = sum;
				}
			}
		}
	}
}

/*********************************************************************
 * _convolveSeparate
 */

static void _convolveSeparate(_KLT_FloatImage imgin,
		ConvolutionKernel horiz_kernel, ConvolutionKernel vert_kernel,
		_KLT_FloatImage imgout)
{
#if 1
	/* Create temporary image */
	_KLT_FloatImage tmpimg;

	tmpimg = _KLTCreateFloatImage(imgin->ncols, imgin->nrows);

	/* Do convolution */
	_convolveImageHoriz(imgin, horiz_kernel, tmpimg);

	_convolveImageVert(tmpimg, vert_kernel, imgout);

	/* Free memory */
	_KLTFreeFloatImage(tmpimg);
#endif
	//_KLTCopyFloatImage(imgin, imgout);
}

/*********************************************************************
 * _convolveSeparate
 */

static void _convolveSeparateB (
  _KLT_FloatImage imgin,
  ConvolutionKernel horiz_kernel,
  ConvolutionKernel vert_kernel,
  _KLT_FloatImage imgout,
  plate pl)
{
  /* Create temporary image */
  _KLT_FloatImage tmpimg;
  tmpimg = _KLTCreateFloatImage(imgin->ncols, imgin->nrows);
  
  //printf("Before _convolveImageHorizB\n");
  /* Do convolution */
  _convolveImageHorizB (imgin, horiz_kernel, tmpimg, pl);

  //printf("Before _convolveImageVertB\n");
  _convolveImageVertB (tmpimg, vert_kernel, imgout, pl);

  /* Free memory */
  _KLTFreeFloatImage(tmpimg);
}

/*********************************************************************
 * _convolveSeparate
 */
static void _convolveSeparateC(vector<KLT_FeatureList> featurelist,
		_KLT_FloatImage imgin, ConvolutionKernel horiz_kernel,
		ConvolutionKernel vert_kernel, _KLT_FloatImage imgout,
		int level, KLT_TrackingContext tc, _KLT_FloatImage tmpimg)
{
	tmpimg->ncols = imgin->ncols;
	tmpimg->nrows = imgin->nrows;
  
	/* Do convolution. */
	_convolveImageHorizC(featurelist, imgin, horiz_kernel, tmpimg,
			level, tc);
	_convolveImageVertC(featurelist, tmpimg, vert_kernel, imgout,
			level, tc);
}

/*********************************************************************
 * _KLTComputeGradients
 */

void _KLTComputeGradients(
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage gradx,
  _KLT_FloatImage grady)
{
				
  /* Output images must be large enough to hold result */
  assert(gradx->ncols >= img->ncols);
  assert(gradx->nrows >= img->nrows);
  assert(grady->ncols >= img->ncols);
  assert(grady->nrows >= img->nrows);

  /* Compute kernels, if necessary */
  if (fabs(sigma - sigma_last) > 0.05)
    _computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
	
  _convolveSeparate(img, gaussderiv_kernel, gauss_kernel, gradx);
  _convolveSeparate(img, gauss_kernel, gaussderiv_kernel, grady);

}
	
/**
 * _KLTComputeGradients
 */
void _KLTComputeGradientsB(_KLT_FloatImage img, float sigma,
		_KLT_FloatImage gradx, _KLT_FloatImage grady, plate pl)
{
	/* Output images must be large enough to hold result. */
	assert(gradx->ncols >= img->ncols);
	assert(gradx->nrows >= img->nrows);
	assert(grady->ncols >= img->ncols);
	assert(grady->nrows >= img->nrows);

	/* Compute kernels, if necessary. */
	if (fabs(sigma - sigma_last) > 0.05) {
		_computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
	}
	_convolveSeparateB(img, gaussderiv_kernel, gauss_kernel, gradx, pl);
	_convolveSeparateB(img, gauss_kernel, gaussderiv_kernel, grady, pl);

}

/**
 * _KLTComputeGradientsC
 */
void _KLTComputeGradientsC(vector<KLT_FeatureList> featurelist,
		_KLT_FloatImage img, float sigma, _KLT_FloatImage gradx,
		_KLT_FloatImage grady, int level, KLT_TrackingContext tc)
{
	_KLT_FloatImage timg;

	/* Output images must be large enough to hold result. */
	assert(gradx->ncols >= img->ncols);
	assert(gradx->nrows >= img->nrows);
	assert(grady->ncols >= img->ncols);
	assert(grady->nrows >= img->nrows);

	timg = _KLTCreateFloatImage(img->ncols, img->nrows);

	/* Compute kernels, if necessary. */
	if (fabs(sigma - sigma_last) > 0.05) {
		_computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
	}

	_convolveSeparateC(featurelist, img, gaussderiv_kernel, gauss_kernel,
			gradx, level, tc, timg);
	_convolveSeparateC(featurelist, img, gauss_kernel, gaussderiv_kernel,
			grady, level, tc, timg);

	_KLTFreeFloatImage(timg);
}

/*********************************************************************
 * _KLTComputeSmoothedImage
 */

void _KLTComputeSmoothedImage(_KLT_FloatImage img, float sigma,
		_KLT_FloatImage smooth)
{
	static int shit = 0;

	/* Output image must be large enough to hold result */
	assert(smooth->ncols >= img->ncols);
	assert(smooth->nrows >= img->nrows);

	/* Compute kernel, if necessary; gauss_deriv is not used */
	if ((fabs(sigma - sigma_last) > 0.05) && (shit == 0)) {
		_computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
		shit = 1;
	}

	_convolveSeparate(img, gauss_kernel, gauss_kernel, smooth);
}

/*********************************************************************
 * _KLTComputeSmoothedImage
 */

void _KLTComputeSmoothedImageB (
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage smooth,
  plate pl)
{
  /* Output image must be large enough to hold result */
  assert(smooth->ncols >= img->ncols);
  assert(smooth->nrows >= img->nrows);

  /* Compute kernel, if necessary; gauss_deriv is not used */
  if (fabs(sigma - sigma_last) > 0.05)
    _computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);

  _convolveSeparateB (img, gauss_kernel, gauss_kernel, smooth, pl);
}

void clearFloatImage(_KLT_FloatImage fimg)
{
	float *ptr;
	int size;

	ptr = fimg->data;
	size = fimg->ncols * fimg->nrows;
	while (size--) {
		*ptr++ = 0.0;
	}
}

/** _KLTComputeSmoothedImage
 */
void _KLTComputeSmoothedImageC(vector<KLT_FeatureList> featurelist,
		_KLT_FloatImage img, float sigma, _KLT_FloatImage smooth,
		int level, KLT_TrackingContext tc)
{
	_KLT_FloatImage timg;

	/* Output image must be large enough to hold result. */
	assert(smooth->ncols >= img->ncols);
	assert(smooth->nrows >= img->nrows);

	timg = _KLTCreateFloatImage(img->ncols, img->nrows);
	clearFloatImage(smooth);

	/* Compute kernel, if necessary; gauss_deriv is not used. */
	if (fabs(sigma - sigma_last) > 0.05) {
		_computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
	}
	_convolveSeparateC(featurelist, img, gauss_kernel, gauss_kernel,
			smooth, level, tc, timg);

	_KLTFreeFloatImage(timg);
}

