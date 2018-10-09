#ifndef __LANCZOS_H__
#define __LANCZOS_H__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "image.h"

double lanczos (int i, int inWidth, int outWidth, double support);

void CalTempContrib (int start, int stop, double *tmpContrib, double *contrib);

int GetRedValue(int rgbValue);

int GetGreenValue(int rgbValue);

int GetBlueValue(int rgbValue);

int ComRGB(int redValue, int greenValue, int blueValue);

int Clip (int x);

int HorizontalFilter (unsigned char *bufImg, int width, int startX, int stopX, int start, int stop, int y, double *pContrib);

int VerticalFilter (unsigned char *pbInImage, int width, int startY, int stopY, int start, int stop, int x, double *pContrib);

unsigned char *HorizontalFiltering (unsigned char *bufImage, int dwInW, int dwInH, int iOutW, int nDots, int nHalfDots, double *contrib, double *tmpContrib,
double *normContrib);

unsigned char *VerticalFiltering (unsigned char *pbImage, int iW, int iH, int iOutH, int nDots, int nHalfDots, double *contrib, double *tmpContrib, double *normContrib);

unsigned char *applylanczos (unsigned char *srcBi, int width, int height, int h, int w_mod, int *W);

#endif
