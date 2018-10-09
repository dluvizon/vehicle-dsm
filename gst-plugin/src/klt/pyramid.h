/**
 * pyramid.h
 */

#ifndef __PYRAMID_H
#define __PYRAMID_H

#include "klt_util.h"

typedef struct {
  int subsampling;
  int nLevels;
  _KLT_FloatImage *img;
  int *ncols, *nrows;
}  _KLT_PyramidRec, *_KLT_Pyramid;


_KLT_Pyramid _KLTCreatePyramid(
  int ncols,
  int nrows,
  int subsampling,
  int nlevels);

_KLT_Pyramid _KLTCopyPyramid (
  _KLT_Pyramid src );

void _KLTComputePyramid(
  _KLT_FloatImage floatimg, 
  _KLT_Pyramid pyramid,
  float sigma_fact);

void _KLTComputePyramidC (
  vector<KLT_FeatureList> featurelist,
  _KLT_FloatImage img,
  _KLT_Pyramid pyramid,
  float sigma_fact,
  KLT_TrackingContext tc);  

void _KLTFreePyramid(
  _KLT_Pyramid pyramid);

#endif /* __PYRAMID_H */

