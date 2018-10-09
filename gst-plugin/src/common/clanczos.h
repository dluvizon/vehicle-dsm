#ifndef _CLANCZOS_H_
#define _CLANCZOS_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "image.h"

double clanczos (int i, int inWidth, int outWidth, double support);

void cCalTempContrib (int start, int stop, double *tmpContrib, double *contrib);

int cGetRedValue(int rgbValue);

int cGetGreenValue(int rgbValue);

int cGetBlueValue(int rgbValue);

int cComRGB(int redValue, int greenValue, int blueValue);

int cClip (int x);

double* cHorizontalFilter (unsigned char *bufImg, int width, int startX, int stopX, int start, int stop, int y, double *pContrib, int nchannels);

double* cVerticalFilter (unsigned char *pbInImage, int width, int startY, int stopY, int start, int stop, int x, double *pContrib, int nchannels);

unsigned char *cHorizontalFiltering (unsigned char *bufImage, int dwInW, int dwInH, int iOutW, int nDots, int nHalfDots, double *contrib, double *tmpContrib,
double *normContrib, int alpha);

unsigned char *cVerticalFiltering (unsigned char *pbImage, int iW, int iH, int iOutH, int nDots, int nHalfDots, double *contrib, double
*tmpContrib, double *normContrib, int alpha);

unsigned char *capplylanczos (unsigned char *srcBi, int width, int height, int h, int w_mod, int *W, int alpha);

#endif
