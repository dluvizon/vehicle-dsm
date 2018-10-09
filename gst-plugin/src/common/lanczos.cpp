#include "lanczos.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double lanczos (int i, int inWidth, int outWidth, double support) {
   double x = (double) i * (double) outWidth / (double) inWidth;
   return sin(x * M_PI) / (x * M_PI) * sin(x * M_PI / support)/ (x * M_PI / support);
}

void CalTempContrib (int start, int stop, double *tmpContrib, double *contrib) {
   double weight = 0;
   int i = 0;
   for (i = start; i <= stop; i++) {
      weight += contrib[i];
   }
   for (i = start; i <= stop; i++) {
      tmpContrib[i] = contrib[i] / weight;
   }
}

int GetRedValue(int rgbValue) {
   return (rgbValue & 0x00ff0000)>>16;
}

int GetGreenValue(int rgbValue) {
   return (rgbValue & 0x0000ff00)>>8;
}

int GetBlueValue(int rgbValue) {
   return rgbValue & 0x000000ff;
}

int ComRGB(int redValue, int greenValue, int blueValue) {
   return (redValue << 16) + (greenValue << 8) + blueValue;
}

int Clip (int x){
   if (x < 0) return 0;
   if (x > 255) return 255;
   return x;
}

int HorizontalFilter (unsigned char *bufImg, int width, int startX, int stopX, int start, int stop, int y, double *pContrib) {
   int valueRGB = 0;
   int i, j;
   for (i = startX, j = start; i <= stopX; i++, j++) {
      valueRGB += bufImg[width * y + i] * pContrib[j];
   }
   return Clip(valueRGB);//ComRGB(Clip((int) valueRed), Clip((int) valueGreen),Clip((int) valueBlue));
}

int VerticalFilter (unsigned char *pbInImage, int width, int startY, int stopY, int start, int stop, int x, double *pContrib) {
   int valueRGB = 0;
   int i, j;
   for (i = startY, j = start; i <= stopY; i++, j++) {
       valueRGB += pbInImage[width * i + x] * pContrib[j];
   }
   return Clip(valueRGB); //ComRGB(Clip((int) valueRed), Clip((int) valueGreen),Clip((int) valueBlue));
   //return valueRGB;
}

unsigned char *HorizontalFiltering (unsigned char *bufImage, int dwInW, int dwInH, int iOutW, int nDots, int nHalfDots, double *contrib, double *tmpContrib, double *normContrib) {
   //int dwInW = bufImage.getWidth();
   //int dwInH = bufImage.getHeight();
   //BufferedImage pbOut = new BufferedImage(iOutW, dwInH, BufferedImage.TYPE_INT_RGB);
   int value = 0;
   unsigned char *pbOut = (unsigned char *)malloc((iOutW * dwInH) * sizeof(unsigned char)); 
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
         CalTempContrib (start, stop, tmpContrib, contrib);
         for (y = 0; y < dwInH; y++) {
            value = HorizontalFilter(bufImage, dwInW, startX, stopX, start, stop, y, tmpContrib);
            //pbOut.setRGB(x, y, value);
            pbOut[y *iOutW  + x] = value;
         }
      } else {
         for (y = 0; y < dwInH; y++) {
            value = HorizontalFilter(bufImage, dwInW, startX, stopX, start, stop, y, normContrib);
            //pbOut.setRGB(x, y, value);
            pbOut[y *iOutW  + x] = value;
         }
      }
   }
   return pbOut;
}

unsigned char *VerticalFiltering (unsigned char *pbImage, int iW, int iH, int iOutH, int nDots, int nHalfDots, double *contrib, double *tmpContrib, double *normContrib) {

   //int iW = pbImage.getWidth();
   //int iH = pbImage.getHeight();
   //BufferedImage pbOut = new BufferedImage(iW, iOutH,BufferedImage.TYPE_INT_RGB);
   int value = 0;
   unsigned char *pbOut = (unsigned char *)malloc((iW * iOutH) * sizeof(unsigned char)); 

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
         CalTempContrib (start, stop, tmpContrib, contrib);
         int x;
         for (x = 0; x < iW; x++) {
            value = VerticalFilter (pbImage, iW, startY, stopY, start, stop, x, tmpContrib);
            //pbOut.setRGB(x, y, value);
            pbOut[y * iW + x] = value;
         }
      } else {
         int x;
         for (x = 0; x < iW; x++) {
            value = VerticalFilter(pbImage, iW, startY, stopY, start, stop, x, normContrib);
            //pbOut.setRGB(x, y, value);
            pbOut[y * iW + x] = value;
         }
      }
   }
   return pbOut;
}

unsigned char *applylanczos (unsigned char *srcBi, int width, int height, int h, int w_mod, int *W) {

  double support = (double) 3.0;
  double scaleV = (double)(h)/(double)(height);
  int w = (int)(floor(width*scaleV/w_mod + 0.5))*w_mod;
  *W = w;
  assert((w % w_mod) == 0);
  double scaleH = (double)(w)/(double)(width);

  if (scaleH >= 1.0 && scaleV >= 1.0) {
     return resize_gray_uchar_image_bilinear (srcBi, width, height, w, h);
  }

  int nHalfDots = (int) ((double) width * support / (double) w);
  int nDots = nHalfDots * 2 + 1;
  double *contrib = (double *)malloc(nDots * sizeof (double)); 
  double *normContrib = (double *)malloc(nDots * sizeof (double));
  double *tmpContrib = (double *)malloc(nDots * sizeof (double));
  int center = nHalfDots;

  if (center < 0) {
     return resize_gray_uchar_image_bilinear (srcBi, width, height, w, h);
  }

  contrib[center] = 1.0;

  double weight = 0.0;
  int i = 0;
  for (i = 1; i <= center; i++) {
     contrib[center + i] = lanczos(i, width,w, support);
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

  unsigned char *pbOut = HorizontalFiltering (srcBi, width, height, w, nDots, nHalfDots, contrib, tmpContrib, normContrib);

  unsigned char *pbFinalOut = VerticalFiltering (pbOut, w, height, h, nDots, nHalfDots, contrib, tmpContrib, normContrib);

  free (pbOut);
  free (contrib);
  free (normContrib);
  free (tmpContrib);

  return pbFinalOut;
}








