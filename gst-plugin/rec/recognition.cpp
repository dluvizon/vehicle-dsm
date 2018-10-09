#include "recognition.h"

void text_recognition (unsigned char *image, int nrows, int ncols, int x, int y, int w, int h, int iframe, int iregion) {

#if 0
  tesseract::TessBaseAPI api;

  api.Init(NULL, "eng");

  api.SetPageSegMode(tesseract::PSM_SINGLE_LINE); /*{PSM_SINGLE_CHAR} is the mode to recognize characters.*/

  int m = 2; /*pixel margin*/

  y -= m;

  x -= m;

  w = (int)(w + 2 * m); /*region width*/

  h = (int)(h + 2 * m); /*region height*/

  PIX *region = pixCreate (w, h, 8);

  /*Cropping the candidate text region from the original image: */
  int i, j, k, l;
  for (j = y, k = 0; j < (y+h); j++, k++) {
     for (i = x, l = 0; i < (x+w); i++, l++) {
        l_uint8 val = image[ncols * j + i]; 
        pixSetPixel (region, l, k, val);
     }
  }

  api.SetImage(region);

  printf("Plate: %s\n", api.GetUTF8Text());

  if (DEBUG_REC) {
     /*Writing the image: */
     char iname[256];
     sprintf(iname, "SAIDA/%05d_%05d.png", iframe, iregion);
     pixWrite (iname, region, IFF_PNG);
     /*Writing the plate: */
     char fname[256];
     sprintf(fname, "SAIDA/%05d_%05d.txt", iframe, iregion);
     FILE *file = fopen (fname, "w");
     fprintf(file, "%s\n", api.GetUTF8Text());
     fclose(file);
  }   
 

  /*tesseract::TessBaseAPI api;

  api.Init(NULL, "eng");

  api.SetPageSegMode(tesseract::PSM_SINGLE_CHAR); */ /*{PSM_SINGLE_CHAR} is the mode to recognize characters.*/ 
#endif
}

