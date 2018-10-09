#include "clanczos.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double clanczos (int i, int inWidth, int outWidth, double support) {
   double PI = (double) 3.14159265358978;
   double x = (double) i * (double) outWidth / (double) inWidth;
   return sin(x * PI) / (x * PI) * sin(x * PI / support)/ (x * PI / support);
}

void cCalTempContrib (int start, int stop, double *tmpContrib, double *contrib) {
   double weight = 0;
   int i = 0;
   for (i = start; i <= stop; i++) {
      weight += contrib[i];
   }
   for (i = start; i <= stop; i++) {
      tmpContrib[i] = contrib[i] / weight;
   }
}

int cGetRedValue(int rgbValue) {
   return (rgbValue & 0x00ff0000)>>16;
}

int cGetGreenValue(int rgbValue) {
   return (rgbValue & 0x0000ff00)>>8;
}

int cGetBlueValue(int rgbValue) {
   return rgbValue & 0x000000ff;
}

int cComRGB(int redValue, int greenValue, int blueValue) {
   return (redValue << 16) + (greenValue << 8) + blueValue;
}

int cClip (int x){
   if (x < 0) return 0;
   if (x > 255) return 255;
   return x;
}

double* cHorizontalFilter (unsigned char *bufImg, int width, int startX, int stopX, int start, int stop, int y, double *pContrib, int nchannels) {

   int i, j, ind;
   
   double *RGB = (double *)malloc(nchannels * sizeof(double));

   for (ind = 0; ind < nchannels; ind++) { RGB[ind] = 0; }

   //printf("maledeto = %d\n", nchannels);
   for (i = startX, j = start; i <= stopX; i++, j++) {
      int ind;
      //printf("%d %d ", y, i);
      for (ind = 0; ind < nchannels; ind++) {
         RGB[ind] += bufImg[nchannels * width * y + nchannels * i + ind] * pContrib[j];
         //printf("%d ", bufImg[nchannels * width * y + nchannels * i + ind]);
         //printf("%f ", RGB[ind]);
      }
      //printf("\n");
   }
   //printf("maledeto 1\n");

   for (ind = 0; ind < nchannels; ind++) { RGB[ind] = cClip(RGB[ind]); }

   return RGB;
}

double* cVerticalFilter (unsigned char *pbInImage, int width, int startY, int stopY, int start, int stop, int x, double *pContrib, int nchannels) {

   int i, j, ind;

   double *RGB = (double *)malloc(nchannels * sizeof(double));

   for (ind = 0; ind < nchannels; ind++) { RGB[ind] = 0; }

   for (i = startY, j = start; i <= stopY; i++, j++) {
      int ind;
      for (ind = 0; ind < nchannels; ind++) {
         RGB[ind] += pbInImage[nchannels * width * i + nchannels * x + ind] * pContrib[j];
      }
   }

   for (ind = 0; ind < nchannels; ind++) { RGB[ind] = cClip(RGB[ind]); }

   return RGB; 
}

unsigned char *cHorizontalFiltering (unsigned char *bufImage, int dwInW, int dwInH, int iOutW, int nDots, int nHalfDots, double *contrib,
double *tmpContrib, double *normContrib, int alpha) {
   
   int nchannels;
   
   if (alpha) { nchannels = 4; }
   
   else { nchannels = 3; }

   unsigned char *pbOut = (unsigned char *)malloc((nchannels * iOutW * dwInH) * sizeof(unsigned char)); 
   int x;
   for (x = 0; x < iOutW; x++) {
      int startX;
      int start;
      int X = (int) (((double) x) * ((double) dwInW) / ((double) iOutW) + 0.5);
      int y = 0;
      startX = X - nHalfDots;
      if (startX < 0) {
         startX = 0;
         start = nHalfDots - X;
      } else {
         start = 0;
      }

      int stop;
      int stopX = X + nHalfDots;
      if (stopX > (dwInW - 1)) {
         stopX = dwInW - 1;
         stop = nHalfDots + (dwInW - 1 - X);
      } else {
         stop = nHalfDots * 2;
      }

      if (start > 0 || stop < nDots - 1) {
         cCalTempContrib (start, stop, tmpContrib, contrib);
         for (y = 0; y < dwInH; y++) {
            //printf("dio mio :D\n");
            double *RGB = cHorizontalFilter(bufImage, dwInW, startX, stopX, start, stop, y, tmpContrib, nchannels);
            //printf("dio mio 1\n");
            int ind;
            for (ind = 0; ind < nchannels; ind++) {
               pbOut[nchannels * y * iOutW + nchannels * x + ind] = (int)(RGB[ind]); 
            }
            //printf("dio mio 2\n");
            free(RGB);
         }
      } else {
         for (y = 0; y < dwInH; y++) {
            //printf("dio mio 00\n");
            double *RGB = cHorizontalFilter(bufImage, dwInW, startX, stopX, start, stop, y, normContrib, nchannels);
            //printf("dio mio 01\n");
            int ind;
            for (ind = 0; ind < nchannels; ind++) {
               pbOut[nchannels * y * iOutW + nchannels * x + ind] = (int)(RGB[ind]); 
            }
            //printf("dio mio 02\n");
            free(RGB);
         }
      }
   }
   return pbOut;
}

unsigned char *cVerticalFiltering (unsigned char *pbImage, int iW, int iH, int iOutH, int nDots, int nHalfDots, double *contrib, double *tmpContrib, double *normContrib, int alpha) {

   int nchannels;
   
   if (alpha) { nchannels = 4; }
   
   else { nchannels = 3; }

   int nch = 3; /*The final image is forced to have 3 channels*/

   unsigned char *pbOut = (unsigned char *)malloc((nch * iW * iOutH) * sizeof(unsigned char)); 

   int y;
   for (y = 0; y < iOutH; y++) {
      int startY;
      int start;
      int Y = (int) (((double) y) * ((double) iH) / ((double) iOutH) + 0.5);
      startY = Y - nHalfDots;
      if (startY < 0) {
         startY = 0;
         start = nHalfDots - Y;
      } else {
         start = 0;
      }
      int stop;
      int stopY = Y + nHalfDots;
      if (stopY > (int) (iH - 1)) {
         stopY = iH - 1;
         stop = nHalfDots + (iH - 1 - Y);
      } else {
         stop = nHalfDots * 2;
      }
      if (start > 0 || stop < nDots - 1) {
         cCalTempContrib (start, stop, tmpContrib, contrib);
         int x;
         for (x = 0; x < iW; x++) {
            double *RGB = cVerticalFilter (pbImage, iW, startY, stopY, start, stop, x, tmpContrib, nchannels);
            int ind;
            for (ind = 0; ind < nch; ind++) {
               pbOut[nch * y * iW + nch * x + ind] = (int)(RGB[ind]); 
            }
            free(RGB);
         }
      } else {
         int x;
         for (x = 0; x < iW; x++) {
            double *RGB = cVerticalFilter (pbImage, iW, startY, stopY, start, stop, x, normContrib, nchannels);
            int ind;
            for (ind = 0; ind < nch; ind++) {
               pbOut[nch * y * iW + nch * x + ind] = (int)(RGB[ind]); 
            }
            free(RGB);
         }
      }
   }
   return pbOut;
}

unsigned char *capplylanczos (unsigned char *srcBi, int width, int height, int h, int w_mod, int *W, int alpha) {

  double support = (double) 3.0;
  double scaleV = (double)(h)/(double)(height);
  int w = (int)(floor(width*scaleV/w_mod + 0.5))*w_mod;
  *W = w;
  assert((w % w_mod) == 0);
  double scaleH = (double)(w)/(double)(width);

  if (scaleH >= 1.0 && scaleV >= 1.0) {
     return resize_color_uchar_image_bilinear (srcBi, width, height, w, h, alpha);
  }

  int nHalfDots = (int) ((double) width * support / (double) w);
  int nDots = nHalfDots * 2 + 1;
  double *contrib = (double *)malloc(nDots * sizeof (double)); 
  double *normContrib = (double *)malloc(nDots * sizeof (double));
  double *tmpContrib = (double *)malloc(nDots * sizeof (double));
  int center = nHalfDots;

  if (center < 0) {
     return resize_color_uchar_image_bilinear (srcBi, width, height, w, h, alpha);
  }

  contrib[center] = 1.0;

  double weight = 0.0;
  int i = 0;
  for (i = 1; i <= center; i++) {
     contrib[center + i] = clanczos(i, width,w, support);
     weight += contrib[center + i];
  }

  for (i = center - 1; i >= 0; i--) {
     contrib[i] = contrib[center * 2 - i];
  }

  weight = weight * 2 + 1.0;

  for (i = 0; i <= center; i++) {
     normContrib[i] = contrib[i] / weight;
  }

  for (i = center + 1; i < nDots; i++) {
     normContrib[i] = normContrib[center * 2 - i];
  }

  unsigned char *pbOut = cHorizontalFiltering (srcBi, width, height, w, nDots, nHalfDots, contrib, tmpContrib, normContrib, alpha);
  unsigned char *pbFinalOut = cVerticalFiltering (pbOut, w, height, h, nDots, nHalfDots, contrib, tmpContrib, normContrib, alpha);

  free (pbOut);
  free (contrib);
  free (normContrib);
  free (tmpContrib);

  return pbFinalOut;
}








