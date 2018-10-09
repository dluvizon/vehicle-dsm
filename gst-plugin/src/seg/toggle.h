#ifndef _TOGGLE_H_
#define _TOGGLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "gtruth/vehicles.h"

#define BACKGROUND 0
#define LOW_VALUE 0
#define HIGH_VALUE 255
#define UNKNOWN_VALUE 0

#define HORZ 98
#define VERT 99

#define LMAX(x,y) (((x) > (y)) ? (x) : (y))
#define LMIN(x,y) (((x) < (y)) ? (x) : (y))

#define PLATE_MIN_WIDTH 18
#define PLATE_MIN_HEIGHT 16

typedef struct plate_ {
   int x;
   int y;
   int w;
   int h;
   float speed;
   struct vehicle_mark *from;
   int lane;
} plate;

void dilate(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *nsrc, unsigned char *ndst,
		unsigned char *isrc, unsigned char *idst,
		int ncols, int nrows, int size, int dir);

void erode(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *nsrc, unsigned char *ndst,
		unsigned char *isrc, unsigned char *idst,
		int ncols, int nrows, int size, int dir);

void min_max (
   int xmin, int ymin, int xmax, int ymax,
   unsigned char *ori_img, unsigned char *ori_min, unsigned char *ori_max,
   unsigned char *inv_img, unsigned char *inv_min, unsigned char *inv_max,
   int ncols, int nrows, int mask_size
);

void toggle(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *ori_img, unsigned char *ori_min,
		unsigned char *ori_max, unsigned char *ori_seg,
		unsigned char *inv_img, unsigned char *inv_min,
		unsigned char *inv_max, unsigned char *inv_seg,
		int ncols, int nrows, int mask_size,
		int contrast, int percentage);

#endif /* _TOGGLE_H_ */
