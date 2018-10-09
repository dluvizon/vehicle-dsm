/**
 * convolve.h
 */

#ifndef __CONVOLVE_H
#define __CONVOLVE_H

#include "klt.h"
#include "klt_util.h"

using namespace std;

void _KLTToFloatImage(const KLT_PixelType *img, int ncols, int nrows,
		_KLT_FloatImage floatimg);

void _KLTComputeGradients (
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage gradx,
  _KLT_FloatImage grady);

void _KLTComputeGradientsB(_KLT_FloatImage img, float sigma,
		_KLT_FloatImage gradx, _KLT_FloatImage grady, plate pl);

void _KLTComputeGradientsC(vector<KLT_FeatureList> featurelist,
		_KLT_FloatImage img, float sigma, _KLT_FloatImage gradx,
		_KLT_FloatImage grady, int level, KLT_TrackingContext tc);

void _KLTGetKernelWidths(float sigma, int *gauss_width,
		int *gaussderiv_width);

void _KLTComputeSmoothedImageC(vector<KLT_FeatureList> featurelist,
		_KLT_FloatImage img, float sigma, _KLT_FloatImage smooth,
		int level, KLT_TrackingContext tc);

void _KLTComputeSmoothedImageB (
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage smooth,
  plate pl);

void _KLTComputeSmoothedImage (
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage smooth);

#endif /* __CONVOLVE_H */

