/**
 * klt_util.h
 */

#ifndef __KLT_UTIL_H
#define __KLT_UTIL_H

typedef struct {
	int ncols;
	int nrows;
	float *data;
}  _KLT_FloatImageRec, *_KLT_FloatImage;

_KLT_FloatImage _KLTCreateFloatImage(
  int ncols, 
  int nrows);

void _KLTCopyFloatImage (
  _KLT_FloatImage src, 
  _KLT_FloatImage dst );

void _KLTFreeFloatImage(
  _KLT_FloatImage);
	
void _KLTPrintSubFloatImage(
  _KLT_FloatImage floatimg,
  int x0, int y0,
  int width, int height);

void _KLTWriteFloatImageToPGM(
  _KLT_FloatImage img,
  const char *filename);

/* for affine mapping */
void _KLTWriteAbsFloatImageToPGM(
  _KLT_FloatImage img,
  const char *filename,float scale);

#endif /* __KLT_UTIL_H */

