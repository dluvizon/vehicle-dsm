/************************************************************************
  Copyright (c) 2003. David G. Lowe, University of British Columbia.
  This software is being made available for research purposes only
  (see file LICENSE for conditions).  This notice must be retained on
  all copies or modifications of this software.
*************************************************************************/

/* sift.c
   This file contains a sample program to extract keypoints from an image
   in PGM format.
*/

#include "key.h"
#include <string.h>

/* -------------------- Local function prototypes ------------------------ */

Image DisplayKeypoints(Image image, Keypoint keys);
void TransformLine(Image image, Keypoint k, float x1, float y1, 
                   float x2, float y2);


/*----------------------------- Routines ----------------------------------*/

/* Top level routine.  Read PGM image from file given in command line 
   arguments, then call GetKeypoints.
*/
int main (int argc, char **argv)
{
    int arg = 0;
    Image image = NULL;
    Keypoint keys = NULL;
    int display = FALSE;


    /* Parse command line arguments and read given files.  The command
       line must specify one input image using command line arguments 
       as follows:
          sift <image.pgm >image.key
          sift -display <image.pgm >result.pgm
    */
    while (++arg < argc) {
      if (! strcmp(argv[arg], "-display")) {
        display = TRUE;
        ++arg;
      }
      else
	FatalError("Invalid command line argument: %s", argv[arg]);
    }

    fprintf(stderr,"Finding keypoints...\n");
    image = ReadPGM(stdin);
    keys = GetKeypoints(image);
    fprintf(stderr,"%d keypoints found.\n", CountKeys(keys));
    if (! display)
      WriteKeys(fileno(stdout),keys);
    else {
      WritePGM(stdout, DisplayKeypoints(image, keys));
      fprintf(stderr,"PGM file output.\n");
    }

    FreeStoragePool(PERM_POOL);
    FreeStoragePool(KEY_POOL);
    exit(0);
}

/* Return a new image that display the keypoints.
*/
Image DisplayKeypoints(Image image, Keypoint keys)
{
    Image result;
    Keypoint k;

    result = CopyImage(image, IMAGE_POOL);

    /* Draw keypoints into result. */
    k = keys;
    while (k != NULL) {
      TransformLine(result, k, 0.0, 0.0, 1.0, 0.0);
      TransformLine(result, k, 0.85, 0.1, 1.0, 0.0);
      TransformLine(result, k, 0.85, -0.1, 1.0, 0.0);
      k = k->next;
    }
    
    return result;
}


/* Draw the given line in the image, but first translate, rotate, and
   scale according to the keypoint parameters.
*/
void TransformLine(Image image, Keypoint k, float x1, float y1, 
                   float x2, float y2)
{
    int rows, cols, len, r1, c1, r2, c2;
    float sine, cosine;

    rows = image->rows;
    cols = image->cols;

    /* The scaling of the unit length arrow is set to approximately the radius
       of the region used to compute the keypoint descriptor. 
    */
    len = 6.0 * k->scale;

    /* Rotate the keypoints by its orientation.
    */
    sine = sin(k->ori);
    cosine = cos(k->ori);

    r1 = (int)(k->row - len * (cosine * y1 + sine * x1) + 0.5);
    c1 = (int)(k->col + len * (- sine * y1 + cosine * x1) + 0.5);
    r2 = (int)(k->row - len * (cosine * y2 + sine * x2) + 0.5);
    c2 = (int)(k->col + len * (- sine * y2 + cosine * x2) + 0.5);

    if (r1 > 0 && r1 < rows && c1 > 0 && c1 < cols &&
        r2 > 0 && r2 < rows && c2 > 0 && c2 < cols)
      DrawLine(image, r1, c1, r2, c2);
}
