#ifndef _RECOGNITION_H_
#define _RECOGNITION_H_


//  #include <tesseract/baseapi.h>
//  #include <leptonica/allheaders.h>
//  #include <tesseract/basedir.h>
//  #include <tesseract/strngs.h>
//  #include <tesseract/tprintf.h>


#define DEBUG_REC 1

  void text_recognition (unsigned char *image, int nrows, int ncols, int x, int y, int w, int h, int iframe, int iregion);

#endif
