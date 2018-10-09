/************************************************************************
  Copyright (c) 2003. David G. Lowe, University of British Columbia.
  This software is being made available for research purposes only
  (see file LICENSE for conditions).  This notice must be retained on
  all copies or modifications of this software.
*************************************************************************/

/* io.c
   This file contains routines for reading and writing PPM files, 
   reading and writing keypoint files, and drawing lines on images.
*/

#include "key.h"
#include <stdarg.h>
#include <zlib.h>


/* -------------------- Local function prototypes ------------------------ */

void SkipComments(FILE *fp);


/*------------------------ Error reporting ----------------------------*/

/* This function prints an error message and exits.  It takes a variable
   number of arguments that function just like those in printf.
*/
void FatalError(char *fmt, ...)
{
    va_list args;

    va_start(args, fmt);
    fprintf(stderr, "Error: ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr,"\n");
    va_end(args);
    exit(1);
}


/* PGM and PPM files allow a comment starting with '#' to end-of-line.  Skip
   white space including any comments.
*/
void SkipComments(FILE *fp)
{
    int ch;

    fscanf(fp," ");      /* Skip white space. */
    while ((ch = fgetc(fp)) == '#') {
      while ((ch = fgetc(fp)) != '\n'  &&  ch != EOF)
	;
      fscanf(fp," ");
    }
    ungetc(ch, fp);      /* Replace last character read. */
}


/*----------------- Read and write PGM files ------------------------*/

/* This reads a PGM file from a given filename and returns the image.
*/
Image ReadPGMFile(char *filename)
{
    Image image;
    FILE *file;

    /* The "b" option is for binary input, which is needed if this is
       compiled under Windows.  It has no effect in Linux.
    */
    file = fopen (filename, "rb");
    if (! file)
	FatalError("Could not open file: %s", filename);

    image = ReadPGM(file);
    fclose(file);

    return image;
}


/* Read a PGM file from the given file pointer and return it as a
   float Image structure with pixels in the range [0,1].   
     See "man pgm" for details on PGM file format.  This handles only
   the usual 8-bit "raw" PGM format.  Use xv or the PNM tools (such as
   pnmdepth) to convert from other formats.
*/
Image ReadPGM(FILE *fp)
{
    int char1, char2, width, height, max, c1, c2, c3, r, c;
    Image image;

    char1 = fgetc(fp);
    char2 = fgetc(fp);
    SkipComments(fp);
    c1 = fscanf(fp, "%d", &width);
    SkipComments(fp);
    c2 = fscanf(fp, "%d", &height);
    SkipComments(fp);
    c3 = fscanf(fp, "%d", &max);

    if (char1 != 'P' || char2 != '5' || c1 != 1 || c2 != 1 || c3 != 1 ||
        max > 255)
      FatalError("Input is not a standard raw 8-bit PGM file.\n"
	    "Use xv or pnmdepth to convert file to 8-bit PGM format.\n");

    fgetc(fp);  /* Discard exactly one byte after header. */

    /* Create floating point image with pixels in range [0,1]. */
    image = CreateImage(height, width, PERM_POOL);
    for (r = 0; r < height; r++)
      for (c = 0; c < width; c++)
        image->pixels[r][c] = ((float) fgetc(fp)) / 255.0;

    return image;
}


/* This writes a PGM file to a given filename in PGM format.
*/
void WritePGMFile(char *filename, Image image)
{
    FILE *file;

    /* The "b" option is for binary input, which is needed if this is
       compiled under Windows.  It has no effect in Linux.
    */
    file = fopen (filename, "wb");
    if (! file)
	FatalError("Could not open file: %s", filename);

    WritePGM(file, image);
    fclose(file);
}


/* Write an image to the file fp in PGM format.
*/
void WritePGM(FILE *fp, Image image)
{
    int r, c, val;

    fprintf(fp, "P5\n%d %d\n255\n", image->cols, image->rows);

    for (r = 0; r < image->rows; r++)
      for (c = 0; c < image->cols; c++) {
	val = (int) (255.0 * image->pixels[r][c]);
	fputc(FUNCMAX(0, FUNCMIN(255, val)), fp);
      }
}


/* Draw a white line from (r1,c1) to (r2,c2) on the image.  Both points
   must lie within the image.
*/
void DrawLine(Image image, int r1, int c1, int r2, int c2)
{
    int i, dr, dc, temp;

    if (r1 == r2 && c1 == c2)  /* Line of zero length. */
      return;

    /* Is line more horizontal than vertical? */
    if (ABS(r2 - r1) < ABS(c2 - c1)) {

      /* Put points in increasing order by column. */
      if (c1 > c2) {
	temp = r1; r1 = r2; r2 = temp;
	temp = c1; c1 = c2; c2 = temp;
      }
      dr = r2 - r1;
      dc = c2 - c1;
      for (i = c1; i <= c2; i++)
	image->pixels[r1 + (i - c1) * dr / dc][i] = 1.0;

    } else {

      if (r1 > r2) {
	temp = r1; r1 = r2; r2 = temp;
	temp = c1; c1 = c2; c2 = temp;
      }
      dr = r2 - r1;
      dc = c2 - c1;
      for (i = r1; i <= r2; i++)
	image->pixels[i][c1 + (i - r1) * dc / dr] = 1.0;
    }
}


/*---------------- Read and write keypoint file ---------------------*/

/* This reads a keypoint file from a given filename and returns the list
   of keypoints.
*/
Keypoint ReadKeyFile(char *filename)
{
    int i, j, num, len;
    Keypoint k, keys = NULL;
    gzFile *fp;

    fp = gzopen (filename, "r");
    if (! fp)
	FatalError("Could not open file: %s", filename);

    if (gzread(fp, &num, sizeof(int)) <= 0)
        FatalError("Invalid keypoint file beginning.");

    if (gzread(fp, &len, sizeof(int)) <= 0)
        FatalError("Invalid keypoint file beginning.");

    if (len != VecLength)
	FatalError("Keypoint descriptor length invalid (should be %d).",
		VecLength);

    for (i = 0; i < num; i++) {
      /* Allocate memory for the keypoint. */
      k = NEW(KeypointSt, KEY_POOL);
      k->next = keys;
      keys = k;
      k->ivec = (unsigned char*) 
     MallocPool(VecLength * sizeof(unsigned char), KEY_POOL);

      if (gzread(fp, &k->row, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      if (gzread(fp, &k->col, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      if (gzread(fp, &k->scale, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      if (gzread(fp, &k->ori, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      for (j = 0; j < len; j++)
	if (gzread(fp, &k->ivec[j], sizeof(char)) <= 0)
	  FatalError("Invalid keypoint file value.");
    }
    gzclose(fp);

    return keys;
}


/* Read keypoints from the given file pointer and return the list of
   keypoints.  The file format starts with 2 integers giving the total
   number of keypoints and the size of descriptor vector for each
   keypoint. Then each keypoint is specified by 4 floating point numbers 
   giving subpixel row and column location, scale, and orientation 
   (in radians from -PI to PI).  Then the descriptor vector for each 
   keypoint is given as a list of integers in range [0,255].

*/
Keypoint ReadKeys(int fd)
{
    int i, j, num, len;
    Keypoint k, keys = NULL;
    gzFile *fp;

    fp = gzdopen (fd, "r");
    if (! fp)
	FatalError("Could not open file stream.");

    if (gzread(fp, &num, sizeof(int)) <= 0)
        FatalError("Invalid keypoint file beginning.");

    if (gzread(fp, &len, sizeof(int)) <= 0)
        FatalError("Invalid keypoint file beginning.");

    if (len != VecLength)
	FatalError("Keypoint descriptor length invalid (should be %d).",
		VecLength);

    for (i = 0; i < num; i++) {
      /* Allocate memory for the keypoint. */
      k = NEW(KeypointSt, KEY_POOL);
      k->next = keys;
      keys = k;
      k->ivec = (unsigned char*) 
     MallocPool(VecLength * sizeof(unsigned char), KEY_POOL);

      if (gzread(fp, &k->row, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      if (gzread(fp, &k->col, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      if (gzread(fp, &k->scale, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      if (gzread(fp, &k->ori, sizeof(float)) <= 0)
        FatalError("Invalid keypoint file format.");

      for (j = 0; j < len; j++)
	if (gzread(fp, &k->ivec[j], sizeof(char)) <= 0)
	  FatalError("Invalid keypoint file value.");
    }
    gzclose(fp);

    return keys;
}


/* Count the number of keys from a given list of keypoints.
*/
int CountKeys(Keypoint keys)
{
    Keypoint k;
    int count;

    k = keys;
    count = 0;
    while (k != NULL) {
      count++;
      k = k->next;
    }

    return(count);
}


/* This writes a list of keypoints to a given filename.
*/
void WriteKeyFile(char *filename, Keypoint keys)
{
    Keypoint k;
    int i;
    gzFile *fp;

    fp = gzopen (filename, "w");
    if (! fp)
	FatalError("Could not open file: %s", filename);

    i = CountKeys(keys);
    gzwrite(fp, &i, sizeof(int));
    i = VecLength;
    gzwrite(fp, &i, sizeof(int));
 
    k = keys;
    while (k != NULL) {
      gzwrite(fp, &k->row, sizeof(float));
      gzwrite(fp, &k->col, sizeof(float));
      gzwrite(fp, &k->scale, sizeof(float));
      gzwrite(fp, &k->ori, sizeof(float));
      for (i = 0; i < VecLength; i++)
        gzwrite(fp, &k->ivec[i], sizeof(char));
      k = k->next;
    }
    gzclose(fp);
}


/* Write a list of keypoints to the file fp.
*/
void WriteKeys(int fd, Keypoint keys)
{
    Keypoint k;
    int i;
    gzFile *fp;

    fp = gzdopen (fd, "w");
    if (! fp)
	FatalError("Could not open file stream.");

    i = CountKeys(keys);
    gzwrite(fp, &i, sizeof(int));
    i = VecLength;
    gzwrite(fp, &i, sizeof(int));
 
    k = keys;
    while (k != NULL) {
      gzwrite(fp, &k->row, sizeof(float));
      gzwrite(fp, &k->col, sizeof(float));
      gzwrite(fp, &k->scale, sizeof(float));
      gzwrite(fp, &k->ori, sizeof(float));
      for (i = 0; i < VecLength; i++)
        gzwrite(fp, &k->ivec[i], sizeof(char));
      k = k->next;
    }
    gzclose(fp);
}
