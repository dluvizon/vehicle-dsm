#include <stdio.h>
#include <stdlib.h>

#include "toggle.h"
#include "utils.h"

void dilate(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *nsrc, unsigned char *ndst,
		unsigned char *isrc, unsigned char *idst,
		int ncols, int nrows, int size, int dir)
{

  int i, j, k;

  unsigned char *nmaxarray = (unsigned char *)malloc((2 * size + 1)*sizeof(unsigned char));
  unsigned char *imaxarray = (unsigned char *)malloc((2 * size + 1)*sizeof(unsigned char));

  int hsize = size / 2;

  if (dir == HORZ) {

     int hsize = size / 2;

     int nsteps = ((xmax - xmin + 1) - 2 * hsize) / size;

     unsigned char *nbuffer = (unsigned char *)malloc((xmax - xmin + 1) * sizeof(unsigned char));

     unsigned char *ibuffer = (unsigned char *)malloc((xmax - xmin + 1) * sizeof(unsigned char));

     for (i = imax (hsize, ymin); i < imin (nrows - hsize, ymax); i++) {

        /*Fill buffer with pixels in byte order: */
        int k = 0;
        for (j = xmin; j < xmax; j++) {
           int p = i * ncols + j;
           nbuffer[k] = nsrc[p]; 
           ibuffer[k] = isrc[p]; 
           k++;
        }

        for (j = 0; j < nsteps; j++) {
           int startmax = (j + 1) * size - 1;
           nmaxarray[size - 1] = nbuffer[startmax];
           imaxarray[size - 1] = ibuffer[startmax];
           for (k = 1; k < size; k++) {
              nmaxarray[size - 1 - k] = LMAX (nmaxarray[size - k + 0], nbuffer[startmax - k]);
              nmaxarray[size - 1 + k] = LMAX (nmaxarray[size + k - 2], nbuffer[startmax + k]);
              imaxarray[size - 1 - k] = LMAX (imaxarray[size - k + 0], ibuffer[startmax - k]);
              imaxarray[size - 1 + k] = LMAX (imaxarray[size + k - 2], ibuffer[startmax + k]);
           }
           int startx = hsize + j * size + xmin;
           ndst[i * ncols + startx] = nmaxarray[0];
           ndst[i * ncols + startx + size - 1] = nmaxarray[2 * size - 2];
           idst[i * ncols + startx] = imaxarray[0];
           idst[i * ncols + startx + size - 1] = imaxarray[2 * size - 2];
           for (k = 1; k < size - 1; k++) {
              int nmaxval = LMAX (nmaxarray[k], nmaxarray[k + size - 1]);
              ndst[i * ncols + startx + k] = nmaxval;
              int imaxval = LMAX (imaxarray[k], imaxarray[k + size - 1]);
              idst[i * ncols + startx + k] = imaxval;
           }
        }
     }
     free (nbuffer);
     free (ibuffer);
  }
  if (dir == VERT) {

     int hsize = size / 2;

     int nsteps = ((ymax - ymin) - 2 * hsize) / size;

     unsigned char *nbuffer = (unsigned char *)malloc((ymax - ymin + 1) * sizeof(unsigned char));
     unsigned char *ibuffer = (unsigned char *)malloc((ymax - ymin + 1) * sizeof(unsigned char));

     for (j = imax (hsize, xmin); j < imin (ncols - hsize, xmax); j++) {

        /*Fill buffer with pixels in byte order: */
        int k = 0;
        for (i = ymin; i < ymax; i++) {
           int p = i * ncols + j;
           nbuffer[k] = nsrc[p]; 
           ibuffer[k] = isrc[p]; 
           k++;
        }

        for (i = 0; i < nsteps; i++) {
           int startmax = (i + 1) * size - 1;
           nmaxarray[size - 1] = nbuffer[startmax];
           imaxarray[size - 1] = ibuffer[startmax];
           for (k = 1; k < size; k++) {
              nmaxarray[size - 1 - k] = LMAX (nmaxarray[size - k], nbuffer[startmax - k]);
              nmaxarray[size - 1 + k] = LMAX (nmaxarray[size + k - 2], nbuffer[startmax + k]);
              imaxarray[size - 1 - k] = LMAX (imaxarray[size - k], ibuffer[startmax - k]);
              imaxarray[size - 1 + k] = LMAX (imaxarray[size + k - 2], ibuffer[startmax + k]);
           }
           int starty = hsize + i * size + ymin;
           ndst[(starty) * ncols + j] = nmaxarray[0];
           ndst[(starty + (size-1)) * ncols + j] = nmaxarray[2 * size - 2];
           idst[(starty) * ncols + j] = imaxarray[0];
           idst[(starty + (size-1)) * ncols + j] = imaxarray[2 * size - 2];
           for (k = 1; k < size - 1; k++) {
              int nmaxval = LMAX (nmaxarray[k], nmaxarray[k + size - 1]);
              ndst[(starty + k) * ncols + j] = nmaxval;
              int imaxval = LMAX (imaxarray[k], imaxarray[k + size - 1]);
              idst[(starty + k) * ncols + j] = imaxval;
           }
        }
     }
     free (nbuffer);
     free (ibuffer);
  }
  free(nmaxarray);
  free(imaxarray);
}

void erode(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *nsrc, unsigned char *ndst,
		unsigned char *isrc, unsigned char *idst,
		int ncols, int nrows, int size, int dir)
{

  int i, j, k;

  unsigned char *nminarray = (unsigned char *)malloc((2 * size + 1)*sizeof(unsigned char));
  unsigned char *iminarray = (unsigned char *)malloc((2 * size + 1)*sizeof(unsigned char));

  int hsize = size / 2;

  if (dir == HORZ) {

     int hsize = size / 2;

     int nsteps = ((xmax - xmin + 1) - 2 * hsize) / size;

     unsigned char *nbuffer = (unsigned char *)malloc((xmax - xmin + 1) * sizeof(unsigned char));

     unsigned char *ibuffer = (unsigned char *)malloc((xmax - xmin + 1) * sizeof(unsigned char));

     for (i = imax (hsize, ymin); i < imin (nrows - hsize, ymax); i++) {

        /*Fill buffer with pixels in byte order: */
        int k = 0;
        for (j = xmin; j < xmax; j++) {
           int p = i * ncols + j;
           nbuffer[k] = nsrc[p]; 
           ibuffer[k] = isrc[p]; 
           k++;
        }

        for (j = 0; j < nsteps; j++) {
           int startmax = (j + 1) * size - 1;
           nminarray[size - 1] = nbuffer[startmax];
           iminarray[size - 1] = ibuffer[startmax];
           for (k = 1; k < size; k++) {
              nminarray[size - 1 - k] = LMIN (nminarray[size - k + 0], nbuffer[startmax - k]);
              nminarray[size - 1 + k] = LMIN (nminarray[size + k - 2], nbuffer[startmax + k]);
              iminarray[size - 1 - k] = LMIN (iminarray[size - k + 0], ibuffer[startmax - k]);
              iminarray[size - 1 + k] = LMIN (iminarray[size + k - 2], ibuffer[startmax + k]);
           }
           int startx = hsize + j * size + xmin;
           ndst[i * ncols + startx] = nminarray[0];
           ndst[i * ncols + startx + size - 1] = nminarray[2 * size - 2];
           idst[i * ncols + startx] = iminarray[0];
           idst[i * ncols + startx + size - 1] = iminarray[2 * size - 2];
           for (k = 1; k < size - 1; k++) {
              int nmaxval = LMIN (nminarray[k], nminarray[k + size - 1]);
              ndst[i * ncols + startx + k] = nmaxval;
              int imaxval = LMIN (iminarray[k], iminarray[k + size - 1]);
              idst[i * ncols + startx + k] = imaxval;
           }
        }
     }
     free (nbuffer);
     free (ibuffer);
  }
  if (dir == VERT) {

     int hsize = size / 2;

     int nsteps = ((ymax - ymin) - 2 * hsize) / size;

     unsigned char *nbuffer = (unsigned char *)malloc((ymax - ymin + 1) * sizeof(unsigned char));
     unsigned char *ibuffer = (unsigned char *)malloc((ymax - ymin + 1) * sizeof(unsigned char));

     for (j = imax (hsize, xmin); j < imin (ncols - hsize, xmax); j++) {

        /*Fill buffer with pixels in byte order: */
        int k = 0;
        for (i = ymin; i < ymax; i++) {
           int p = i * ncols + j;
           nbuffer[k] = nsrc[p]; 
           ibuffer[k] = isrc[p]; 
           k++;
        }

        for (i = 0; i < nsteps; i++) {
           int startmax = (i + 1) * size - 1;
           nminarray[size - 1] = nbuffer[startmax];
           iminarray[size - 1] = ibuffer[startmax];
           for (k = 1; k < size; k++) {
              nminarray[size - 1 - k] = LMIN (nminarray[size - k], nbuffer[startmax - k]);
              nminarray[size - 1 + k] = LMIN (nminarray[size + k - 2], nbuffer[startmax + k]);
              iminarray[size - 1 - k] = LMIN (iminarray[size - k], ibuffer[startmax - k]);
              iminarray[size - 1 + k] = LMIN (iminarray[size + k - 2], ibuffer[startmax + k]);
           }
           int starty = hsize + i * size + ymin;
           ndst[(starty) * ncols + j] = nminarray[0];
           ndst[(starty + (size-1)) * ncols + j] = nminarray[2 * size - 2];
           idst[(starty) * ncols + j] = iminarray[0];
           idst[(starty + (size-1)) * ncols + j] = iminarray[2 * size - 2];
           for (k = 1; k < size - 1; k++) {
              int nmaxval = LMIN (nminarray[k], nminarray[k + size - 1]);
              ndst[(starty + k) * ncols + j] = nmaxval;
              int imaxval = LMIN (iminarray[k], iminarray[k + size - 1]);
              idst[(starty + k) * ncols + j] = imaxval;
           }
        }
     }
     free (nbuffer);
     free (ibuffer);
  }
  free(nminarray);
  free(iminarray);
}

void min_max (
     int xmin, int ymin, int xmax, int ymax, 
     unsigned char *ori_img, unsigned char *ori_min, unsigned char *ori_max, 
     unsigned char *inv_img, unsigned char *inv_min, unsigned char *inv_max, 
     int ncols, int nrows, int mask_size
  ) {

   int i, j, k, l;

   int hmask = mask_size / 2;

   int nsup, ninf, isup, iinf;
  
   for (i = imax (hmask, ymin); i < imin (nrows - hmask, ymax); i++) {
      for (j = imax (hmask, xmin); j < imin (ncols - hmask, xmax); j++) {
         int pi = i * ncols + j;
         nsup = ninf = ori_img[pi];
         isup = iinf = inv_img[pi];
         for (k = -hmask; k <= hmask; k++) {
            for (l = -hmask; l <= hmask; l++) {
               int pm = (i + k) * ncols + (j + l);
               if (nsup < ori_img[pm]) {
                  nsup = ori_img[pm];
               }
               else if (ninf > ori_img[pm]) {
                  ninf = ori_img[pm];
               }
               if (isup < inv_img[pm]) {
                  isup = inv_img[pm];
               }
               else if (iinf > inv_img[pm]) {
                  iinf = inv_img[pm];
               }
            }
         }
         ori_min[pi] = ninf; ori_max[pi] = nsup;
         inv_min[pi] = iinf; inv_max[pi] = isup;
      }
   }
}

void toggle(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *ori_img, unsigned char *ori_min,
		unsigned char *ori_max, unsigned char *ori_seg,
		unsigned char *inv_img, unsigned char *inv_min,
		unsigned char *inv_max, unsigned char *inv_seg,
		int ncols, int nrows, int mask_size,
		int contrast, int percentage)
{
	int x, y;
	int p;
	int hmask = mask_size / 2;

	for (y = imax(hmask, ymin); y < imin(nrows - hmask, ymax); y++) {
		for (x = imax(hmask, xmin);
				x < imin(ncols - hmask, xmax); x++) {
			p = y * ncols + x;

			/* Segmentation at the original image: */
			if ((ori_max[p] - ori_min[p]) < contrast) {
				ori_seg[p] = UNKNOWN_VALUE;
			} else if ((ori_max[p] - ori_min[p]) < contrast) {
				ori_seg[p] = UNKNOWN_VALUE;
			} else {
				if ((ori_max[p] - ori_img[p]) <	percentage *
						(ori_img[p] -ori_min[p])/100) {
					ori_seg[p] = HIGH_VALUE;
				} else {
					ori_seg[p] = LOW_VALUE;
				}
			}

			/* Segmentation at the inverse image: */
			if ((inv_max[p] - inv_min[p]) < contrast ) {
				inv_seg[p] = UNKNOWN_VALUE;
			} else if ((inv_max[p] - inv_min[p]) < contrast) {
				inv_seg[p] = UNKNOWN_VALUE;
			} else {
				if ((inv_max[p] - inv_img[p]) < percentage *
						(inv_img[p] -inv_min[p])/100) {
					inv_seg[p] = HIGH_VALUE;
				} else {
					inv_seg[p] = LOW_VALUE;
				}
			}
		}
	}
}

