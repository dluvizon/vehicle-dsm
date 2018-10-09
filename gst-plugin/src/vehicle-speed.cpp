/**
 * @file vehicle-speed.c
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 29/01/2014
 */

#define _DLEVEL VDSM_DLEVEL

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <locale.h>

#include "vehicle-speed.h"
#include "dsm-common.h"
#include "index.h"
#include "font.h"
#include "video-defs.h"

#include "klt/klt.h"
#include "klt/pyramid.h"
#include "sift/key.h"
#include "seg/toggle.h"
#include "snoopertext/snoopertext.h"
#include "snoopertext/detection-minetto.h"
#include "snoopertext/detection-zheng.h"
//#include "swt/swt-wrapper.h"
#include "svm/svm.h"
#include "speed/speed.h"
#include "thog/thog.h"
#include "debug.h"

#define WINDOW_LEFT_MARGIN 24
#define WINDOW_RIGHT_MARGIN 24
#define WINDOW_TOP_MARGIN 280
#define WINDOW_BOTTON_MARGIN 24

#define MINY_TRACKING 220

#define WINDOW_SHIFT_X_MAX 12
#define WINDOW_SHIFT_Y_MAX 12

#if (WINDOW_SHIFT_X_MAX >= WINDOW_LEFT_MARGIN)
# error WINDOW_LEFT_MARGIN must be greater than WINDOW_SHIFT_X_MAX
#endif
#if (WINDOW_SHIFT_X_MAX >= WINDOW_RIGHT_MARGIN)
# error WINDOW_RIGHT_MARGIN must be greater than WINDOW_SHIFT_X_MAX
#endif
#if (WINDOW_SHIFT_Y_MAX >= WINDOW_BOTTON_MARGIN)
# error WINDOW_BOTTON_MARGIN must be greater than WINDOW_RIGHT_MARGIN
#endif

const bool use_detection = true;
static char *old_locale;

volatile int _dbg_tab;

static int total_slopes;
static int total_slopes_gt;
static int slopes_motocycle;
static int slopes_motocycle_gt;
static int slopes_car;
static int slopes_car_gt;

struct show_speed {
	int cnt_show;
	int show_m;
	float speed;
	float measured;
	float sum_meas;
	int cnt_meas;
	bool wait_meas;
	int x;
	int y;
};

# define NUM_SHOW_SPEED 3
static struct show_speed sspeed[NUM_SHOW_SPEED];

/*------------------------ Error reporting ----------------------------*/
/* This function prints an error message and exits.  It takes a variable
   number of arguments that function just like those in printf.
*/
void FatalError(const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	fprintf(stderr, "Error: ");
	vfprintf(stderr, fmt, args);
	fprintf(stderr,"\n");
	va_end(args);
	exit(1);
}

static int dsm_select_lane(plate p, struct lane_limits *lanes)
{
	int x;
	int y;

	x = p.x + p.w / 2;
	y = p.y + p.h / 2;
	if (x < 0) {
		x = 0;
	}
	if (x >= lanes->width) {
		x = lanes->width - 1;
	}
	if (y < 0) {
		y = 0;
	}
	if (y >= lanes->height) {
		y = lanes->height - 1;
	}

	return lanes->lmask[x + y * lanes->width];
}

static void dsm_select_plates(vector<plate> &plates,
		struct lane_limits *lanes, int iframe)
{
	int plate_maxy[GTRUTH_MAX_LANES] = {0};
	int plate_k[GTRUTH_MAX_LANES] = {0};
	vector<plate> tmp;
	int max_w;
	int max_h;
	int lane;
	int x;
	int y;
	int k;

	max_w = lanes->width - 2;
	max_h = lanes->height - 2;

	/* Remove invalid plates or plates too small. */
	for (k = 0; k < plates.size();) {
		if (plates.at(k).w < PLATE_MIN_WIDTH) {
			plates.at(k).x -=
				(PLATE_MIN_WIDTH - plates.at(k).w) / 2;
			plates.at(k).w = PLATE_MIN_WIDTH;
		}
		if (plates.at(k).h < PLATE_MIN_HEIGHT) {
			plates.at(k).y -=
				(PLATE_MIN_HEIGHT - plates.at(k).h) / 2;
			plates.at(k).h = PLATE_MIN_HEIGHT;
		}
		if ((plates.at(k).x < 0) ||
				(plates.at(k).y < 0) ||
				(plates.at(k).x + plates.at(k).w > max_w) ||
				(plates.at(k).y + plates.at(k).h > max_h)) {
			plates.erase(plates.begin() + k);
			continue;
		}
		k++;
	}

	/* Remove plates outside the lane's boudary. */
	for (k = 0; k < plates.size();) {
		lane = dsm_select_lane(plates.at(k), lanes);

		/*
		printf("          p.x / w = %d / %d\n",
				plates.at(k).x, plates.at(k).w);
		printf("          p.y / h = %d / %d\n",
				plates.at(k).y, plates.at(k).h);
		printf("          lane    = %d\n", lane);
		*/

		if ((lane < 1) || (lane > GTRUTH_MAX_LANES)) {
			/* This plate is out of the lane's boundary. */
			plates.erase(plates.begin() + k);
		} else {
			plates.at(k).lane = lane;
			k++;
		}
	}
	/* Find the highest plate in Y axis indise a single lane. */
	for (k = 0; k < plates.size(); k++) {
		x = plates.at(k).x + plates.at(k).w / 2;
		y = plates.at(k).y + plates.at(k).h / 2;
		if (x < 0) {
			x = 0;
		}
		if (x >= lanes->width) {
			x = lanes->width - 1;
		}
		if (y < 0) {
			y = 0;
		}
		if (y >= lanes->height) {
			y = lanes->height - 1;
		}
		/* Swip lane from 1~MAX to 0~(MAX-1). */
		lane = lanes->lmask[x + y * lanes->width] - 1;
		if (y > plate_maxy[lane]) {
			plate_maxy[lane] = y;
			plate_k[lane] = k;
		}
	}
	/* Copy the highest plates in each lane. */
	for (k = 0; k < GTRUTH_MAX_LANES; k++) {
		if (0 == plate_maxy[k]) {
			continue;
		}
		tmp.push_back(plates.at(plate_k[k]));
	}
	plates.clear();

	/* Copy the selected plates back to the plates vector. */
	for (k = 0; k < tmp.size(); k++) {
		/*
		printf("          p.x = %d\n", tmp.at(k).x);
		printf("          p.y = %d\n", tmp.at(k).y);
		printf("          p.w = %d\n", tmp.at(k).w);
		printf("          p.h = %d\n", tmp.at(k).h);
		printf("          lanes w,h = %d,%d\n",
				lanes->width, lanes->height);
		*/
		tmp.at(k).speed = 0.0;
		tmp.at(k).from = NULL;
		plates.push_back(tmp.at(k));
	}
}

static void dsm_paint_plates(uint8_t *img, int w, int h, vector<plate> &plates,
		uint8_t cb, uint8_t cr)
{
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	int k;

	for (k = 0; k < plates.size(); k++) {
		xmin = plates.at(k).x;
		xmax = xmin + plates.at(k).w;
		ymin = plates.at(k).y;
		ymax = ymin + plates.at(k).h;

		if (ymin < 1) {
			ymin = 1;
		}
		if (ymax > h - 2) {
			ymax = h - 2;
		}
		if (xmin < 1) {
			xmin = 1;
		}
		if (xmax > w - 2) {
			xmax = w - 2;
		}
		paint_chrom_rect(img, w, h, xmin, ymin, xmax, ymax, cb, cr);
		paint_lum_rect(img, w, h, xmin - 2, ymin - 2,
				xmax + 2, ymax + 2, 5);

		/*
		draw_text(img, w, xmin - 42, ymin - 78, 4, 255, 0, " LICENSE ");
		draw_text(img, w, xmin - 42, ymin - 44, 4, 255, 0, "  PLATE  ");
		paint_chrom_rect(img, w, h, xmin - 42, ymin - 78,
				xmin + 246, ymin - 7, 72, 220);
				*/
	}

}

static void check_lane_already_tracking(struct lane_limits *lanes,
		vector<KLT_FeatureList> &features, int miny)
{
	unsigned int i;
	int lane;
	int x, y;

	/* Clear tracking history. */
	for (i = 0; i < GTRUTH_MAX_LANES; i++) {
		lanes->area[i].tracking_plate = false;
	}

	/**
	 * For each feature, calculates the license plate center point and
	 * check in what lane this license plate is in.  After, if the Y value
	 * is greather than 'miny', than set the flag tracking_plate to true.
	 */
	for (i = 0; i < features.size(); i++) {
		x = (features.at(i)->xmax + features.at(i)->xmin) / 2;
		y = (features.at(i)->ymax + features.at(i)->ymin) / 2;
		assert((x >= 0) && (x < lanes->width) &&
				(y >= 0) && (y < lanes->height));
		lane = lanes->lmask[x + y * lanes->width];
		if ((lane < 1) || (lane > GTRUTH_MAX_LANES)) {
			continue;
		}
		lane--;
		if (y > miny) {
			lanes->area[lane].tracking_plate = true;
		}
	}
}

static void check_slope_already_tracking(vector<KLT_FeatureList> &features,
		vector<struct sub_slope> &slopes)
{
	int k;
	int i;

	for (k = 0; k < slopes.size(); k++) {
		int r1x1 = slopes.at(k).left;
		int r1x2 = slopes.at(k).right;
		int r1y1 = slopes.at(k).up;
		int r1y2 = slopes.at(k).down;

		for (i = 0; i < features.size(); i++) {
			int r2x1 = features.at(i)->xmin;
			int r2x2 = features.at(i)->xmax;
			int r2y1 = features.at(i)->ymin;
			int r2y2 = features.at(i)->ymax;

			/* Testing overlap: */
			if ((r1x1 < r2x2) && (r2x1 < r1x2) &&
					(r1y1 < r2y2) && (r2y1 < r1y2)) {
				slopes.at(k).search_plate = false;
			}
		}
	}
}

#if 0
void remove_duplicates(vector<KLT_FeatureList> &features,
		vector<struct sub_slope> &slopes, int iframe)
{
	int i;

   for (i = 0; i < slopes.size(); i++) {
      int slope_xmin = slopes.at(i).left;
      int slope_ymin = slopes.at(i).up;
      int slope_xmax = slopes.at(i).right;
      int slope_ymax = slopes.at(i).down;
      int k;
      int no_plates = 0, plate_pos = 0, plate_ymax = 0;
      for (k = 0; k < features.size(); k++) {
         int j;
         int xmin = 99999, ymin = 99999, xmax = 0, ymax = 0;
         int in = 0;
         for (j = 0; (j < features.at(k)->nFeatures); j++) {
             if (features.at(k)->feature[j]->val >= 0 && (features.at(k)->discovered + 3 < iframe) ) {
                if (features.at(k)->feature[j]->x < xmin) { xmin = features.at(k)->feature[j]->x; }  
                if (features.at(k)->feature[j]->y < ymin) { ymin = features.at(k)->feature[j]->y; }  
                if (features.at(k)->feature[j]->x > xmax) { xmax = features.at(k)->feature[j]->x; }  
                if (features.at(k)->feature[j]->y > ymax) { ymax = features.at(k)->feature[j]->y; }  
                in = 1;
             }
         }
         if ( (xmin >= slope_xmin) && (ymin >= slope_ymin) && (xmax <= slope_xmax) && (ymax <= slope_ymax) && in) {
            if (no_plates == 0) {
               plate_pos = k;
               plate_ymax = ymax;
            } 
            else { 
               if (plate_ymax < ymax) {
                  features.at(plate_pos)->falsepositive = 1;
                  plate_ymax = ymax;
                  plate_pos = k;  
               }
               else {
                  features.at(k)->falsepositive = 1;
               } 
            }  
            no_plates++;
         } 
      }
   }
}
#endif

static void vs_unnecessary_remove_features(vector<KLT_FeatureList> &features)
{
	vector<KLT_FeatureList>::iterator it;
	KLT_FeatureList element;
	int i, j, k;

	/* for all features list, do ... */
	for (it = features.begin(); it != features.end();) {
		int in = 0;

		element = *it;
		for (j = 0; j < element->nFeatures; j++) {
			if (element->feature[j]->val >= 0) {
				in++;
			}
		}

                if (element->falsepositive || (in == 0)) {
			KLTFreeFeatureList(element);
			it = features.erase(it);
		} else {
			it++;
		}
	}
}

extern "C" {

Image crop_image(const unsigned char *src, int nrows, int ncols,
		int xmin, int ymin, int xmax, int ymax)
{
   if (xmin < 0) { xmin = 0; }
   if (ymin < 0) { ymin = 0; }
   if (xmax >= ncols) { xmax = ncols - 1; }
   if (ymax >= nrows) { ymax = nrows - 1; }

   Image dst = CreateImage (ymax - ymin + 1, xmax - xmin + 1, PERM_POOL);

   for (int r = 0; r < dst->rows; r++) {
      for (int c = 0; c < dst->cols; c++) {
         int p = ncols * (r + ymin) + (c + xmin);
         dst->pixels[r][c] = ((float)src[p]) / 255.0;
      }
   }
   return dst;
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

/* Return a new image that contains the two images with im1 above im2.
*/
Image CombineImagesVertically(Image im1, Image im2)
{
    int rows, cols, r, c;
    Image result;

    rows = im1->rows + im2->rows;
    cols = FUNCMAX(im1->cols, im2->cols);
    result = CreateImage(rows, cols, PERM_POOL);

    /* Set all pixels to 0,5, so that blank regions are grey. */
    for (r = 0; r < rows; r++)
      for (c = 0; c < cols; c++)
        result->pixels[r][c] = 0.5;

    /* Copy images into result. */
    for (r = 0; r < im1->rows; r++)
      for (c = 0; c < im1->cols; c++)
        result->pixels[r][c] = im1->pixels[r][c];
    for (r = 0; r < im2->rows; r++)
      for (c = 0; c < im2->cols; c++)
        result->pixels[r + im1->rows][c] = im2->pixels[r][c];

    return result;
}

/* Return squared distance between two keypoint descriptors.
*/
int DistSquared (Keypoint k1, Keypoint k2)
{
    int i, dif, distsq = 0;
    unsigned char *pk1, *pk2;

    pk1 = k1->ivec;
    pk2 = k2->ivec;

    for (i = 0; i < 128; i++) {
      dif = (int) *pk1++ - (int) *pk2++;
      distsq += dif * dif;
    }
    return distsq;
}

/* This searches through the keypoints in klist for the two closest
   matches to key.  If the closest is less than 0.6 times distance to
   second closest, then return the closest match.  Otherwise, return
   NULL.
*/
Keypoint CheckForMatch(Keypoint key, Keypoint klist)
{
    int dsq, distsq1 = 100000000, distsq2 = 100000000;
    Keypoint k, minkey = NULL;

    /* Find the two closest matches, and put their squared distances in
       distsq1 and distsq2.
    */
    for (k = klist; k != NULL; k = k->next) {
      dsq = DistSquared(key, k);

      if (dsq < distsq1) {
        distsq2 = distsq1;
        distsq1 = dsq;
        minkey = k;
      } else if (dsq < distsq2) {
        distsq2 = dsq;
      }
    }

    /* Check whether closest distance is less than 0.6 of second. */
    //if (10 * 10 * distsq1 < 8 * 8 * distsq2)
    if (10 * 10 * distsq1 < 6 * 6 * distsq2)
      return minkey;
    else return NULL;
}

static void WritePGM(FILE *fp, Image image)
{
    int r, c, val;

    fprintf(fp, "P5\n%d %d\n255\n", image->cols, image->rows);

    for (r = 0; r < image->rows; r++)
      for (c = 0; c < image->cols; c++) {
	val = (int) (255.0 * image->pixels[r][c]);
	fputc(FUNCMAX(0, FUNCMIN(255, val)), fp);
      }
}

static void WritePGMFile(char *filename, Image image)
{
    FILE *file;

    /* The "b" option is for binary input, which is needed if this is
       compiled under Windows.  It has no effect in Linux.
    */
    file = fopen(filename, "wb");
    if (! file)
	FatalError("Could not open file: %s", filename);

    WritePGM(file, image);
    fclose(file);
}

/**
 * Given a pair of images and their keypoints, pick the first keypoint from one
 * image and find its closest match in the second set of keypoints. Then write
 * the result to a file.
 */
extern "C" void FindMatches(Image im1, Keypoint keys1, Image im2,
		Keypoint keys2, KLT_FeatureList fl,
		int shift_x1, int shift_y1, int shift_x2, int shift_y2,
		int iframe, const char *output_path)
{
	Keypoint k;
	Keypoint match;
	Image result1;
	Image result2;
	int count = 0;

	/* Create a new image that joins the two images vertically. */
        #if VDSM_DLEVEL > 3
	   result1 = CombineImagesVertically(im2, im1);
	   result2 = CombineImagesVertically(im2, im1);
        #endif

	/* Match the keys in list keys1 to their best matches in keys2. */
	fl->dx = 0.0;
	fl->dy = 0.0;

        /* TESTE PARA REMOCAO DE OUTLIER */
        int outliers = 0;
        static const int SDEV = 2;
        static const double CDEV = 0.5;
        do {
           outliers = 0;
           int nfeatures = 0;
           double sum_x = 0.0, sum_y = 0.0;
           for (k= keys1; k != NULL; k = k->next) {
              match = CheckForMatch(k, keys2);
              if ((match != NULL) && (k->valid)) {
                 double d_x = ((fl->xmin - shift_x2) + match->col) - ((fl->xmin - shift_x1) + k->col);
                 double d_y = ((fl->ymin - shift_y2) + match->row) - ((fl->ymin - shift_y1) + k->row);
                 sum_x += d_x;
                 sum_y += d_y;
                 nfeatures++;
              } 
           }
           double mu_x = (sum_x / (double)nfeatures);
           double mu_y = (sum_y / (double)nfeatures);
           //printf("Mu - %f %f\n", mu_x, mu_y);
           sum_x = 0.0; sum_y = 0.0;
           for (k= keys1; k != NULL; k = k->next) {
              match = CheckForMatch(k, keys2);
              if ((match != NULL) && (k->valid)) {
                 double d_x = ((fl->xmin - shift_x2) + match->col) - ((fl->xmin - shift_x1) + k->col);
                 double d_y = ((fl->ymin - shift_y2) + match->row) - ((fl->ymin - shift_y1) + k->row);
                 sum_x += (d_x - mu_x)*(d_x - mu_x);
                 sum_y += (d_y - mu_y)*(d_y - mu_y);
              } 
           }
           double dev_x = sqrt(sum_x / (double) (nfeatures - 1));
           double dev_y = sqrt(sum_y / (double) (nfeatures - 1));
           //printf("Dev - %f %f\n", dev_x, dev_y);
           for (k= keys1; k != NULL; k = k->next) {
              match = CheckForMatch(k, keys2);
              if ((match != NULL) && (k->valid)) {
                 double d_x = ((fl->xmin - shift_x2) + match->col) - ((fl->xmin - shift_x1) + k->col);
                 double d_y = ((fl->ymin - shift_y2) + match->row) - ((fl->ymin - shift_y1) + k->row);
                 if ( (fabs(d_x - mu_x) > SDEV * dev_x && (dev_x > CDEV)) || ((fabs(d_y - mu_y) > SDEV * dev_y) && (dev_y > CDEV)) ) {
                    k->valid = 0;
                    outliers++;
                 }
                 /*else {
                    printf("Inlier - %f %f\n", dev_x, dev_y);
   		    fl->dx += ((fl->xmin - shift_x2) + match->col) - ((fl->xmin - shift_x1) + k->col);
	            fl->dy += ((fl->ymin - shift_y2) + match->row) - ((fl->ymin - shift_y1) + k->row);
     	            count++;
                    if (DEBUG_SIFT) {
                        printf("Desenho - %f %f\n", dev_x, dev_y);
		        DrawLine(result, (int) k->row, (int) k->col,(int) (match->row + im1->rows),(int) match->col);
                    }
                 }*/
              } 
           }
        } while (outliers > 0); 

	for (k= keys1; k != NULL; k = k->next) {
		match = CheckForMatch(k, keys2);

		/* Draw all keypoints matching. */
		if (match != NULL) {
                        #if VDSM_DLEVEL > 3
			    DrawLine(result1, (int) match->row, (int) match->col,
			    		     (int) (k->row + im2->rows),
					     (int) k->col);
                        #endif
		}
		/* Draw a line on the image from keys1 to match.  Note that we
		 * must add row count of first image to row position in second
		 * so that line ends at correct location in second image.
		 */
		if ((match != NULL) && (k->valid)) {
			count++;
                        #if VDSM_DLEVEL > 3
			    DrawLine(result2, (int) match->row, (int) match->col,
			    		     (int) (k->row + im2->rows),
					     (int) k->col);
                        #endif
			/*
			fl->dy += (match->row - (180 + k->row));
			fl->dx += (match->col - (180 + k->col));
			fl-dy += (match->row - k->row);
			*/
			fl->dy += ((fl->ymin - shift_y2) + match->row) -
				((fl->ymin - shift_y1) + k->row);
			fl->dx += ((fl->xmin - shift_x2) + match->col) -
				((fl->xmin - shift_x1) + k->col);
		}
	}
	if (count != 0) {
		fl->dx /= count;
		fl->dy /= count;
	} else {
		fl->dx = 0.0;
		fl->dy = 0.0;
	}

        #if VDSM_DLEVEL > 3
	   {
		   /* Write the SIFT matches is enabled. */
		   char _fname[256];
		   sprintf(_fname, "/tmp/sift/f%05d_r1.pgm", iframe);
		   WritePGMFile(_fname, result1);
		   sprintf(_fname, "/tmp/sift/f%05d_r2.pgm", iframe);
		   WritePGMFile(_fname, result2);
	   }
	   DisallocMatrix(result1->pixels, result1->rows, result1->cols, PERM_POOL);
	   DisallocMatrix(result2->pixels, result2->rows, result2->cols, PERM_POOL);
	   free(result1);
	   free(result2);
        #endif
}

void update_displacement (KLT_FeatureList fl) {

      int j;

      double mean_velx = 0.0;
      double mean_vely = 0.0;

      int nFeat = 0;

      fl->xmin = fl->ymin = +9999999;
      fl->xmax = fl->ymax = -9999999;

      for (j = 0; j < fl->nFeatures; j++) {

         if (fl->feature[j]->val == KLT_TRACKED) {

            mean_velx += fl->feature[j]->x - fl->feature[j]->prev_x;
            mean_vely += fl->feature[j]->y - fl->feature[j]->prev_y;

            if (fl->feature[j]->x < fl->xmin) {
               fl->xmin = fl->feature[j]->x;
            }
            if (fl->feature[j]->x > fl->xmax) {
               fl->xmax = fl->feature[j]->x;
            }
            if (fl->feature[j]->y < fl->ymin) {
               fl->ymin = fl->feature[j]->y;
            }
            if (fl->feature[j]->y > fl->ymax) {
               fl->ymax = fl->feature[j]->y;
            }

            nFeat++;
         }
      }
      fl->dx = mean_velx/nFeat;
      fl->dy = mean_vely/nFeat;
}

static void draw_circle_nv12(uint8_t *chrom, int w, int h, int x, int y,
		int cb, int cr)
{
	print_pixel_nv12(chrom, w, h, x, y, cb, cr);
	print_pixel_nv12(chrom, w, h, x + 2, y, cb, cr);
	print_pixel_nv12(chrom, w, h, x - 2, y, cb, cr);
	print_pixel_nv12(chrom, w, h, x, y + 2, cb, cr);
	print_pixel_nv12(chrom, w, h, x, y - 2, cb, cr);
}

void draw_circle_lum(uint8_t *lum, int w, int h, int x, int y,
		int gray)
{
	print_pixel_lum(lum, w, h, x, y, 255 - gray);

	print_pixel_lum(lum, w, h, x + 1, y, 255 - gray);
	print_pixel_lum(lum, w, h, x - 1, y, 255 - gray);
	print_pixel_lum(lum, w, h, x, y + 1, 255 - gray);
	print_pixel_lum(lum, w, h, x, y - 1, 255 - gray);

	print_pixel_lum(lum, w, h, x + 2, y, gray);
	print_pixel_lum(lum, w, h, x - 2, y, gray);
	print_pixel_lum(lum, w, h, x, y + 2, gray);
	print_pixel_lum(lum, w, h, x, y - 2, gray);

	print_pixel_lum(lum, w, h, x + 1, y + 1, gray);
	print_pixel_lum(lum, w, h, x + 1, y - 1, gray);
	print_pixel_lum(lum, w, h, x - 1, y + 1, gray);
	print_pixel_lum(lum, w, h, x - 1, y - 1, gray);
}

void draw_line_lum(uint8_t *lum, int y1, int x1, int y2, int x2,
		int w, int h, int gray, int thickness)
{
	int temp;
	int dr;
	int dc;
	int i;
	int t;

	/* line of zero length */
	if (y1 == y2 && x1 == x2) {
		return;
	}
	/* is line more horizontal than vertical? */
	if (ABS(y2 - y1) < ABS(x2 - x1)) {

		/* put points in increasing order by column */
		if (x1 > x2) {
			temp = y1;
			y1 = y2;
			y2 = temp;
			temp = x1;
			x1 = x2;
			x2 = temp;
		}
		dr = y2 - y1;
		dc = x2 - x1;
		for (i = x1; i <= x2; i++) {
			int y = y1 + (i - x1) * dr / dc;
			int x = i;

			for (t = y - thickness / 2; t < y + (thickness + 1) / 2; t++) {
				print_pixel_lum(lum, w, h, x, t, gray);
			}
		}
	} else {
		if (y1 > y2) {
			temp = y1;
			y1 = y2;
			y2 = temp;
			temp = x1;
			x1 = x2;
			x2 = temp;
		}
		dr = y2 - y1;
		dc = x2 - x1;
		for (i = y1; i < y2; i++) {
			int y = i;
			int x = x1 + (i - y1) * dc / dr;

			for (t = x - thickness / 2; t < x + (thickness + 1) / 2; t++) {
				print_pixel_lum(lum, w, h, t, y, gray);
			}
		}
	}
}

static void draw_line_nv12(uint8_t *chrom, int y1, int x1, int y2, int x2,
		int w, int h, int cb, int cr)
{
	int temp;
	int dr;
	int dc;
	int i;

	/* line of zero length */
	if (y1 == y2 && x1 == x2) {
		return;
	}
	/* is line more horizontal than vertical? */
	if (ABS(y2 - y1) < ABS(x2 - x1)) {

		/* put points in increasing order by column */
		if (x1 > x2) {
			temp = y1;
			y1 = y2;
			y2 = temp;
			temp = x1;
			x1 = x2;
			x2 = temp;
		}
		dr = y2 - y1;
		dc = x2 - x1;
		for (i = x1; i <= x2; i++) {
			int y = y1 + (i - x1) * dr / dc;
			int x = i;

			print_pixel_nv12(chrom, w, h, x, y, cb, cr);
		}
	} else {
		if (y1 > y2) {
			temp = y1;
			y1 = y2;
			y2 = temp;
			temp = x1;
			x1 = x2;
			x2 = temp;
		}
		dr = y2 - y1;
		dc = x2 - x1;
		for (i = y1; i < y2; i++) {
			int y = i;
			int x = x1 + (i - y1) * dc / dr;

			print_pixel_nv12(chrom, w, h, x, y, cb, cr);
		}
	}
}

void draw_motion_vectors_lum(uint8_t *image, int width, int height,
		int iframe, KLT_FeatureList features)
{
	int min_x = INT_MAX;
	int min_y = INT_MAX;
	int act_x;
	int act_y;
	int prv_x;
	int prv_y;

	for (int i = 0; i < features->nFeatures; i++)  {
		if (features->feature[i]->val == KLT_TRACKED) {
			act_x = (int) features->feature[i]->prv_x[(iframe) % SMPSIZE];
			act_y = (int) features->feature[i]->prv_y[(iframe) % SMPSIZE];
			prv_x = (int) features->feature[i]->prv_x[(iframe - 1) % SMPSIZE];
			prv_y = (int) features->feature[i]->prv_y[(iframe - 1) % SMPSIZE];
			if ((act_x != NOT_FILLED) && (act_y != NOT_FILLED) &&
					(prv_x != NOT_FILLED) && (prv_y != NOT_FILLED)) {
				draw_line_lum(image, prv_y, prv_x, act_y,
						act_x, width, height, 255, 1);
			}
			act_x = (int) features->feature[i]->x;
			act_y = (int) features->feature[i]->y;
			if (act_x < min_x) {
				min_x = act_x;
			}
			if (act_y < min_y) {
				min_y = act_y;
			}
			draw_circle_lum(image, width, height, act_x, act_y, 255);
		}
	}
}

static void draw_motion_vectors_across(uint8_t *image, int width, int height,
		int frame_i, KLT_FeatureList features, int thickness)
{
	int act_x, act_y, prv_x, prv_y;
	int MAX_INT = 99999;
	int min_x = MAX_INT, min_y = MAX_INT;
	uint8_t *chrom;

	chrom = image + width * height;
	for (int i = 0; i < features->nFeatures; i++)  {
		if (features->feature[i]->val == KLT_TRACKED) {
			for (int k = 0; k < (SMPSIZE - 1); k++) {
				act_x = (int) features->feature[i]->
					prv_x[(frame_i - k) % SMPSIZE];
				act_y = (int) features->feature[i]->
					prv_y[(frame_i - k) % SMPSIZE];
				prv_x = (int) features->feature[i]->
					prv_x[(frame_i - k - 1) % SMPSIZE];
				prv_y = (int) features->feature[i]->
					prv_y[(frame_i - k - 1) % SMPSIZE];
				if ((act_x != NOT_FILLED) &&
						(act_y != NOT_FILLED) &&
						(prv_x != NOT_FILLED) &&
						(prv_y != NOT_FILLED) ) {
					draw_line_lum(image, prv_y, prv_x,
							act_y, act_x, width,
							height, 244 - 28 * k,
							thickness);
					draw_line_nv12(chrom, prv_y, prv_x, act_y,
						act_x, width, height,
						240, 64 + 16 * k);
						//64, 255 - 36 * k);
				}
			}
			act_x = (int) features->feature[i]->x;
			act_y = (int) features->feature[i]->y;
			if (act_x < min_x) {
				min_x = act_x;
			}
			if (act_y < min_y) {
				min_y = act_y;
			}
			//draw_circle_lum(image, width, height, act_x, act_y, 255);
		}
		if (features->feature[i]->val == KLT_TRACKED) {
			for (int k = 0; k < (SMPSIZE - 1); k++) {
				act_x = (int) features->feature[i]->
					prv_x[(frame_i - k) % SMPSIZE];
				act_y = (int) features->feature[i]->
					prv_y[(frame_i - k) % SMPSIZE];
				if ((act_x != NOT_FILLED) &&
						(act_y != NOT_FILLED)) {
					draw_circle_lum(image, width, height, act_x, act_y, 255);
				}
			}
		}
	}
}

static void draw_motion_vectors_nv12(uint8_t *image, int width, int height,
		int frame_i, KLT_FeatureList features, int color)
{
	uint8_t *chrom;
	int act_x, act_y, prv_x, prv_y;
	int MAX_INT = 99999;
	int min_x = MAX_INT, min_y = MAX_INT;

	chrom = image + width * height;
	for (int i = 0; i < features->nFeatures; i++)  {
		if (features->feature[i]->val == KLT_TRACKED) {
			act_x = (int) features->feature[i]->x;
			act_y = (int) features->feature[i]->y;
			if (act_x < min_x) {
				min_x = act_x;
			}
			if (act_y < min_y) {
				min_y = act_y;
			}
			draw_circle_nv12(chrom, width, height, act_x, act_y,
					color, 255 - color);

			for (int k = 0; k < (SMPSIZE - 1); k++) {
				act_x = (int) features->feature[i]->
					prv_x[(frame_i - k) % SMPSIZE];
				act_y = (int) features->feature[i]->
					prv_y[(frame_i - k) % SMPSIZE];
				prv_x = (int) features->feature[i]->
					prv_x[(frame_i - k - 1) % SMPSIZE];
				prv_y = (int) features->feature[i]->
					prv_y[(frame_i - k - 1) % SMPSIZE];
				if ((act_x != NOT_FILLED) &&
						(act_y != NOT_FILLED) &&
						(prv_x != NOT_FILLED) &&
						(prv_y != NOT_FILLED) ) {
					draw_line_nv12(chrom, prv_y, prv_x, act_y,
						act_x, width, height,
						color, 255 - 42 * k);
				}
			}
		}
	}
}

extern int CountKeys(Keypoint keys);

struct sub_data *initialize_background_module(const char *opath, int w, int h)
{
	char strtmp[256];
	struct sub_data *sub;
	struct rectangle_m win = {
		.xmin = WINDOW_LEFT_MARGIN,
		.ymin = WINDOW_TOP_MARGIN,
		.xmax = w - WINDOW_RIGHT_MARGIN,
		.ymax = h - WINDOW_BOTTON_MARGIN,
	};
	size_t len;

	sub = background_sub_new(w, h, win,
			SUBSAMPLE_STEP_X, SUBSAMPLE_STEP_Y);
	if (!sub) {
		exit(-1);
	}
	sub->iframe_last = 0;

	len = strlen(opath);
	if (len > sizeof(strtmp)) {
		print_err("in _%s_ string output_path is too large", __func__);
	}
	strncpy(strtmp, opath, sizeof(strtmp));
	strncat(strtmp, "/backgroundsub", sizeof(strtmp) - len);

	return sub;
}

/**
 * Read and parse the XML given by lanes_xml and fill the structs area,
 * build the polygons poly and measure.
 */
int vs_set_lanes_area(struct lane_area *area,
		const char *lanes_xml)
{
	pugi::xml_document doc;
	pugi::xml_parse_result result;
	pugi::xml_node root;
	pugi::xml_node lanes;
	pugi::xml_node nLanes;
	pugi::xml_node it;
	size_t nlanes;
	int ret = 0;

	/* Initialize the xml lib with the given file. */
	result = doc.load_file(lanes_xml);
	root = doc.child("RegionOfInterest");
	lanes = root.child("lanes");
	nLanes = root.child("nLanes");

	nlanes = atoi(nLanes.attribute("total").value());
	if (nlanes > GTRUTH_MAX_LANES) {
		print_err("Too many lanes in file %s! Max is %d.\n",
				lanes_xml, GTRUTH_MAX_LANES);
		return -1;
	}
	memset(area, 0, nlanes * sizeof(struct lane_area));

	for (it = lanes.first_child(); it; it = it.next_sibling()) {
		struct lane_area *lp;
		struct polygon_4s *poly;
		pugi::xml_node p;
		int id;

		id = atoi(it.attribute("id").value());
		if ((id < 1) || (id > nlanes)) {
			return -1;
		}
		lp = area + (id - 1);

		if (strstr(it.name(), "lane")) {
			poly = &lp->poly;
		} else if (strstr(it.name(), "measure")) {
			poly = &lp->measure;
		} else {
			/* Error, skip it. */
			print_warn("Element not recognized %s\n", it.name());
			continue;
		}

		p = it.child("pTopLeft");
		poly->p[0].x = atof(p.attribute("x").value());
		poly->p[0].y = atof(p.attribute("y").value());
		p = it.child("pTopRight");
		poly->p[1].x = atof(p.attribute("x").value());
		poly->p[1].y = atof(p.attribute("y").value());
		p = it.child("pBotLeft");
		poly->p[2].x = atof(p.attribute("x").value());
		poly->p[2].y = atof(p.attribute("y").value());
		p = it.child("pBotRight");
		poly->p[3].x = atof(p.attribute("x").value());
		poly->p[3].y = atof(p.attribute("y").value());

		printf("%f %f\n", poly->p[0].x, poly->p[0].y);
		printf("%f %f\n", poly->p[1].x, poly->p[1].y);
		printf("%f %f\n", poly->p[2].x, poly->p[2].y);
		printf("%f %f\n\n", poly->p[3].x, poly->p[3].y);

		ret = compute_lines_polygon_4s(poly);
		if (-1 == ret) {
			return -1;
		}
	}

	return 0;
}

/** Vehicle speed preload function.
 * This function must be called at the first
 * frame. The image width and height are assumed to be constant to all the
 * next frames.
 * @param buff Pointer to the first image frame.
 * @param width Image width.
 * @param height Image height.
 * @return A pointer to a newly allocated structure dsm_vs_data.
 */
void *vs_preload(struct dsm_common_data *cd, int width, int height)
{
	struct dsm_vs_data *vs;
	char fname[256];
	int nr_class;
	int size;
	int ret;
	int i;

	old_locale = strdup(setlocale(LC_ALL, NULL));
	setlocale(LC_ALL, "C");

	print_dbg("Vehicle Detection Speed is preloading...\n");
	vs = (struct dsm_vs_data *) malloc(sizeof(struct dsm_vs_data));
	if (!vs) {
		goto end;
	}
	vs->img_width = width;
	vs->img_height = height;
	vs->object_counter = 0;
	cd->timestamp = 0.0;

	snprintf(fname, sizeof(fname), "%s/../../matrix.txt", cd->output_path);
	ret = load_ipm_matrix_from_file(fname);
	if (ret) {
		exit(-1);
	}

	vs->background_sub = initialize_background_module(cd->output_path,
			width, height);

	vs->vt = gtruth_vehicles_new(cd->vehicles_xml);

	snprintf(fname, sizeof(fname), "%s/../../lanes.xml", cd->output_path);
	ret = vs_set_lanes_area(vs->lanes.area, fname);
	if (-1 == ret) {
		exit(-1);
	}
	vs->lanes.lmask = (uint8_t *)
		calloc(1, width * height * sizeof(uint8_t));
	vs->lanes.mmask = (uint8_t *)
		calloc(1, width * height * sizeof(uint8_t));
	vs->lanes.width = width;
	vs->lanes.height = height;
	for (i = 0; i < GTRUTH_MAX_LANES; i++) {
		build_mask_from_polygon(vs->lanes.lmask, width, height,
				&vs->lanes.area[i].poly, i + 1);
		build_mask_from_polygon(vs->lanes.mmask, width, height,
				&vs->lanes.area[i].measure, i + 1);
		vs->lanes.area[i].tracking_plate = false;
	}

	vs->tc = KLTCreateTrackingContext();
	vs->tc_gt = KLTCreateTrackingContext();
	vs->last = 0;
	vs->kltbufs.tmp = _KLTCreateFloatImage(width, height);
	vs->kltbufs.aux = _KLTCreateFloatImage(width, height);

	print_dbg("loading svm model\n");
	ret = snoopertext_load_svm_model(&vs->snooper, SVM_MODEL_FILE);
	if (-1 == ret) {
		print_err("could not load file " SVM_MODEL_FILE);
		exit(-1);
	}
	print_dbg("loading thog settings\n");
	ret = snoopertext_load_thog_settings(&vs->snooper, THOG_SETTINGS_FILE);
	if (-1 == ret) {
		print_err("could not load file " THOG_SETTINGS_FILE);
		exit(-1);
	}
	print_dbg("alocating memory for snoopertext\n");
	size = width * height;
	vs->snooper.nmin = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.nmax = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.imin = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.imax = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.nseg = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.iseg = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.ntmp = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.itmp = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.invert = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.edges = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.enhanced = (unsigned char *) malloc(size * sizeof(char));
	vs->snooper.dilate = (unsigned char *) malloc(size * sizeof(char));
	nr_class = svm_get_nr_class(vs->snooper.model);
	vs->snooper.prob = (double *) malloc(nr_class * sizeof(double));
	vs->should_save_frame = false;

	total_slopes = 0;
	total_slopes_gt = 0;
	slopes_motocycle = 0;
	slopes_motocycle_gt = 0;
	slopes_car = 0;
	slopes_car_gt = 0;

	for (i = 0; i < NUM_SHOW_SPEED; i++) {
		sspeed[i].cnt_show = 0;
		sspeed[i].show_m = 0;
		sspeed[i].y = 890;
		sspeed[i].wait_meas = false;
	}
	sspeed[0].x = 60;
	sspeed[1].x = 730;
	sspeed[2].x = 1500;

	print_dbg("vs_preload done\n");

end:
	return (void *) vs;
}

#if (VDSM_DLEVEL > 1)
static void vs_write_features_to_image(uint8_t *buff, int w, int h,
		KLT_FeatureList fl)
{
	uint8_t *chrom;
	int x, y;
	int mi;
	int f;

	chrom = buff + w * h;
	for (f = 0; f < fl->nFeatures; f++) {
		if (fl->feature[f]->val >= 0) {
			x = (int) fl->feature[f]->x;
			y = (int) fl->feature[f]->y;
#if (VDSM_DLEVEL > 2)
			mi = INDEX_INIT(w, x, y);
			buff[mi] = 0;
			buff[mi + 1] = 0;
			buff[mi - 1] = 0;
			buff[mi + w] = 0;
			buff[mi - w] = 0;
			buff[mi + w + 1] = 255;
			buff[mi + w - 1] = 255;
			buff[mi - w + 1] = 255;
			buff[mi - w - 1] = 255;
			buff[mi + 2] = 255;
			buff[mi - 2] = 255;
			buff[mi + 2 * w] = 255;
			buff[mi - 2 * w] = 255;
#endif
			print_pixel_nv12(chrom, w, h, x, y, 255, 255);
			print_pixel_nv12(chrom, w, h, x + 2, y + 2, 255, 255);
			print_pixel_nv12(chrom, w, h, x + 2, y - 2, 255, 255);
			print_pixel_nv12(chrom, w, h, x - 2, y + 2, 255, 255);
			print_pixel_nv12(chrom, w, h, x - 2, y - 2, 255, 255);
		}
	}
}
#else
#define vs_write_features_to_image(...)
#endif

static void vs_print_vehicle_plates(uint8_t *frame, struct dsm_vs_data *vs,
		long int iframe, int thickness)
{
	struct vehicle_table *vt;
	struct vehicle_mark *mark;
	uint8_t *chrom;
	int width;
	int height;
	int i, j;

	vt = vs->vt;
	if (iframe >= vt->size) {
		/* Out of range. */
		return;
	}
	if (vt->vehicles[iframe]) {
		/* If exist a vehicle in this frame. */
		mark = vt->vehicles[iframe];
		width = vs->img_width;
		height = vs->img_height;
		chrom = frame + (width * height);

		for (j = mark->y; j < mark->y + mark->h; j += 2) {
			for (i = mark->x; i < mark->x + mark->w; i += 2) {
				print_pixel_nv12(chrom, width, height,
						i, j, 0, 255);
			}
		}
		paint_lum_rect(frame, width, height, mark->x - 2, mark->y - 2,
				mark->x + mark->w + 2, mark->y + mark->h + 2,
				thickness);
	}
}

static void vs_print_lanes_to_nv12(uint8_t *chrom, int w, int h, uint8_t *mask,
		int sum)
{
	uint8_t cb;
	uint8_t cr;
	int x, y;
	int pix;

	for (y = 0; y < h; y += 2) {
		for (x = 0; x < w; x += 2) {
			pix = mask[y * w + x];
			if (0 == pix) {
				continue;
			}
			if (1 == pix) {
				cb = 122 + sum;
				cr = 64 - sum;
			} else if (2 == pix) {
				cb = 104 - sum;
				cr = 128 + sum;
			} else if (3 == pix) {
				cb = 142 + sum;
				cr = 156 + sum;
			}
			print_pixel_nv12(chrom, w, h, x, y, cb, cr);
		}
	}
}

void vs_log_vehicle_speed(KLT_FeatureList fl, int iframe,
		const char *out_path, const char *log_type, bool *save_frame)
{
	char fname[256];
	FILE *file;
	double s;
	double r;
	int xav;
	int yav;
	int ret;
	float sp, er;
	int remain_feat;

	sp = fl->speed[iframe % SMPSIZE];
	if (fl->vehicle.radar) {
		er = sp - fl->vehicle.speed;
	} else {
		er = 0.0;
	}

	yav = (fl->ymin + fl->ymax) / 2;
	xav = (fl->xmin + fl->xmax) / 2;

	remain_feat = KLTCountRemainingFeatures(fl);
	if (0 == remain_feat) {
		/* Nothing to do if there isn't features tracked. */
		return;
	}

	snprintf(fname, sizeof(fname), "%s/%s/f%d/%04ld.txt",
			out_path, log_type, fl->vehicle.lane, fl->object_id);
	ret = access(fname, F_OK);
	if (-1 == ret) {
		/* the file does not exist */
		file = fopen(fname, "w");
	} else {
		file = fopen(fname, "a");
	}
	if (NULL == file) {
		printf("Error: could not open file '%s'\n", fname);
		exit(1);
	}
	if (-1 == ret) {
		/* if the files was created now, put the header first */
		fprintf(file, "# Log file generated by " PACKAGE_STRING "\n\n");
		fprintf(file, "# Global parameters:\n");
		fprintf(file, "id : %d\n", fl->object_id);
		fprintf(file, "lane : %d\n", fl->vehicle.lane);
		fprintf(file, "plate : %d\n", !!fl->vehicle.plate);
		fprintf(file, "radar : %d\n", !!fl->vehicle.radar);
		fprintf(file, "sema : %d\n", !!fl->vehicle.sema);
		fprintf(file, "moto : %d\n", !!fl->vehicle.moto);
		fprintf(file, "real speed : %.1f\n", fl->vehicle.speed);
		fprintf(file, "l.p. width : %d\n", fl->vehicle.w);
		fprintf(file, "l.p. height : %d\n", fl->vehicle.h);
		fprintf(file, "# ---------------------------------------\n\n");

		fprintf(file,   "frame\t"
				"speed\t"
				"err\t"
				"x\t"
				"y\t"
				"dx\t"
				"dy\t"
				"ipm_x\t"
				"ipm_y\t"
				"nfeat\n");
		fprintf(file, "# ---------------------------------------------"
				"--------------------------------\n");
	}
	fprintf(file, "%05d\t", iframe);
	fprintf(file, "%.1f\t", sp);
	fprintf(file, "%.1f\t", sp - fl->vehicle.speed);
	fprintf(file, "%d\t", xav);
	fprintf(file, "%d\t", yav);
	fprintf(file, "%.1f\t", fl->dx);
	fprintf(file, "%.1f\t", fl->dy);
	fprintf(file, "%.1f\t", fl->ipm_x);
	fprintf(file, "%.1f\t", fl->ipm_y);
	fprintf(file, "%d\n", remain_feat);
	fclose(file);

	if (fl->vehicle.radar && (fl->vehicle.speed > 10.0) &&
			(sp > 10.0) && (fabs(er) < 10.0)) {
		if ((er > 2.0) | (er < -3.0)) {
			*save_frame = true;
		}
	}
	if (!(*save_frame) && (remain_feat < 2)) {
		*save_frame = true;
	}
}

static void search_groundtruth_from_xml(vector<plate> &plates,
		struct vehicle_table *vt, int iframe)
{
	const int search_frame_before = 15;
	const int search_frame_after = 10;
	int frame_max;
	int frame_min;
	int j, k;

	if ((iframe >= vt->size) || (iframe < 1)) {
		return;
	}
	frame_min = max(0, iframe - search_frame_before);
	frame_max = min(vt->size - 1, iframe + search_frame_after);

	for (k = 0; k < plates.size(); k++) {
		for (j = frame_min; j <= frame_max; j++) {
			if (!vt->vehicles[j]) {
				continue;
			}
			if ((vt->vehicles[j]->lane != plates.at(k).lane)) {
				continue;
			}
			if (vt->vehicles[j]->used) {
				continue;
			}
			vt->vehicles[j]->used = true;
			plates.at(k).from = vt->vehicles[j];
		}
	}
}

static void save_frame_window_to_pgm(uint8_t *img, int imgw, int w, int h,
		const char *filename)
{
	FILE *file;
	int pix_cnt;
	int mi;
	int x;
	int y;

	file = fopen(filename, "w");
	if (!file) {
		print_err("could not open file %s", filename);
		return;
	}

	pix_cnt = 0;
	mi = 0;

	fprintf(file, "P2\n");
	fprintf(file, "%d %d\n", w, h);
	fprintf(file, "255\n");

	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			fprintf(file, "%d ", (int) img[mi]);
			if (pix_cnt % 12 == 0) {
				fprintf(file, "\n");
			}
			pix_cnt++;
			mi++;
		}
		mi += (imgw - w);
	}

	fclose(file);
}

static void save_slopes_to_pgm(uint8_t *img, int w, int h, const char *path,
		vector<struct sub_slope> &slopes, int iframe)
{
	char filename[256];
	int x1, x2;
	int y1, y2;
	int ret;
	int k;

	for (k = 0; k < slopes.size(); k++) {
		x1 = slopes.at(k).left;
		x2 = slopes.at(k).right;
		y1 = slopes.at(k).up;
		y2 = slopes.at(k).down;
		if ((x1 < 0) || (x1 >= x2) || (x2 >= w) ||
				(y1 < 0) || (y1 >= y2) || (y2 >= h)) {
			fprintf(stderr, "Error: save_slopes: "
					"slope out of range\n");
			continue;
		}
		ret = snprintf(filename, sizeof(filename),
				"%s/slopes/%05d_s%02d.pgm", path, iframe, k);
		if (ret >= sizeof(filename)) {
			fprintf(stderr, "Error: save_slopes: "
					"filename too large\n");
			return;
		}
		save_frame_window_to_pgm(img + x1 + (y1 * w), w,
				x2 - x1, y2 - y1, filename);
	}
}

static vector<plate> read_plates_from_xml(struct vehicle_table *vt, int iframe)
{
	vector<plate> plates;
	struct vehicle_mark *vehicle;
	plate p;

	if ((iframe >= vt->size) || (iframe < 0)) {
		/* Out of range for plates. */
		goto end;
	}
	if (!vt->vehicles[iframe]) {
		/* If there is not a vehicle in this position, return. */
		goto end;
	}
	vehicle = vt->vehicles[iframe];
	p.x = vehicle->x;
	p.y = vehicle->y;
	p.w = vehicle->w;
	p.h = vehicle->h;
	if (vehicle->radar) {
		p.speed = vehicle->speed;
	} else {
		p.speed = 0.0;
	}
	p.from = vehicle;
	plates.push_back(p);

end:
	return plates;
}

static int lane_by_plate(plate p)
{
	int x;

	x = p.x + p.w / 2;
	if (x < 534) {
		return 1;
	}
	if (x < 1300) {
		return 2;
	}
	return 3;
}


void vs_end_of_stream_handler(void *_data)
{
	print_warn("Exiting software...\n");

	/* Configure the old locale back. */
	setlocale(LC_ALL, old_locale);
	free(old_locale);

	exit(0);
}


static void draw_speed_values(uint8_t *img, int width, int height,
		int iframe, KLT_FeatureList fl, const char *out_path)
{
	int x, y;

	x = (fl->xmin + fl->xmax) / 2 - 100;
	y = (fl->ymin + fl->ymax) / 2 - 50;
	if ((x < 10) || (x > width - 100)) {
		return;
	}
	if ((y < 40) || (y > height - 40)) {
		return;
	}
	if (fl->vehicle.radar) {
		draw_text(img, width, x, y, 3, 0, 255, " %.1f / %.1f",
				fl->speed[iframe % SMPSIZE], fl->vehicle.speed);
	} else {
		draw_text(img, width, x, y, 3, 0, 255, " %.1f",
				fl->speed[iframe % SMPSIZE]);
	}

        write_nv12_to_jpeg(img, width, height, "%s/snapshot/vel_%05d.jpg",
			out_path, iframe);
}

static void decolorize(uint8_t *img, int width, int height)
{
	uint8_t *chr;
	int size;

	chr = img + width * height;
	size = width * height / 2;
	while (size--) {
		*chr = (*chr + 128) / 2;
		chr++;
	}
}

KLT_FeatureList match_plate_with_feature_list(plate *pplate,
		vector<KLT_FeatureList> &features)
{
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	int k;

	for (k = 0; k < features.size(); k++) {
		KLT_FeatureList fl = features.at(k);
		xmin = max(fl->xmin, pplate->x);
		xmax = min(fl->xmax, pplate->x + pplate->w);
		ymin = max(fl->ymin, pplate->y);
		ymax = min(fl->ymax, pplate->y + pplate->h);
		if ((xmin < xmax) && (ymin < ymax)) {
			/* These rectangles are matching. */
			return fl;
		}
	}

	return NULL;
}

/**
 * Para cada placa no array plates, verifica se o retângulo dessa placa
 * coincide com algum conjunto de features sendo rastreados.  Se sim, apenas
 * copia as informações da placa para o conjunto de features.  Se não, cria um
 * novo conjunto de features na região da placa usando o KLT e atribui a esse
 * conjunto as informações contidas na placa.
 */
void compute_features_from_plates(vector<plate> &plates,
		vector<KLT_FeatureList> &features, struct dsm_vs_data *vs,
		struct dsm_common_data *cd, bool atrack)
{
	uint8_t *roimg;
	int width;
	int height;
	int k;

	if (0 == plates.size()) {
		return;
	}
	width = vs->img_width;
	height = vs->img_height;
	roimg = grab_past_buffer(cd, 0);

	for (k = 0; k < plates.size(); k++) {
		KLT_FeatureList fl;
		/**
		 * First of all, check if this plate is in the same region as
		 * some feature list.  If it is, just copy the vehicle
		 * structure from this plate to the corresponding feature list.
		 */
		fl = match_plate_with_feature_list(&plates.at(k), features);
		if (!fl) {
			/* If not, create a new feature list. */
			fl = KLTCreateFeatureList(10);
			KLTSelectGoodFeatures(vs->tc, roimg, width, height,
					fl, plates.at(k));

			fl->discovered = cd->iframe;
			fl->falsepositive = 0;
			fl->object_id = ++vs->object_counter;
			fl->print_speed = false;
			fl->print_measured = false;
			/*
			printf("  DBG: nova featureList criada: frame %d id %d\n",
					fl->discovered, fl->object_id);
					*/
			features.push_back(fl);
		}
		if (plates.at(k).from) {
			fl->vehicle = *plates.at(k).from;
		} else {
			/* Null vehicle pointer from. */
			fl->vehicle.iframe = cd->iframe;
			fl->vehicle.lane =
				dsm_select_lane(plates.at(k), &vs->lanes);
			fl->vehicle.plate = true;
			fl->vehicle.sema = false;
			fl->vehicle.moto = false;
			fl->vehicle.radar = false;
			fl->vehicle.x = plates.at(k).x;
			fl->vehicle.y = plates.at(k).y;
			fl->vehicle.w = plates.at(k).w;
			fl->vehicle.h = plates.at(k).h;
			fl->vehicle.frame_start = cd->iframe;
			fl->vehicle.frame_end = cd->iframe;
			fl->vehicle.speed = 0.0;
		}
		if (fl->vehicle.plate) {
			fl->always_track = atrack;
		} else {
			fl->always_track = false;
		}
	}
}

/**
 * Return the percentage of mathin of the recond rectangle r2 againts
 * the first rectangle r1.
 * The return value is from 0.0. to 1.0.
 */
float match_rectangles(int r1x1, int r1x2, int r1y1, int r1y2,
		int r2x1, int r2x2, int r2y1, int r2y2)
{
	float a_match;
	float a_ref;
	int x1;
	int x2;
	int y1;
	int y2;

	if (r1x2 < r1x1 || r1y2 < r1y1 || r2x2 < r2x1 || r2y2 < r2y1) {
		printf("\n  Skipping due to some shit...\n\n");
		return 0.0;
	}

	x1 = max(r1x1, r2x1);
	x2 = min(r1x2, r2x2);
	y1 = max(r1y1, r2y1);
	y2 = min(r1y2, r2y2);
	if (x1 >= x2) {
		return 0.0;
	}
	if (y1 >= y2) {
		return 0.0;
	}
	a_match = (x2 - x1 + 1) * (y2 - y1 + 1);
	a_ref = (r1x2 - r1x1 + 1) * (r1y2 - r1y1 + 1);

	return a_match / a_ref;
}

struct vehicle_mark *find_vehicle_mark_in_slope(struct sub_slope *s,
		vector<KLT_FeatureList> &f)
{
	float match;
	int k;

	for (k = 0; k < f.size(); k++) {
		/**
		 * Match the slope region with each feature list bounding box.
		 */
		match = match_rectangles(f.at(k)->xmin, f.at(k)->xmax,
				f.at(k)->ymin, f.at(k)->ymax,
				s->left, s->right, s->up, s->down);
		if (match > 0.5) {
			return &f.at(k)->vehicle;
		}
	}

	return NULL;
}

void debug_guess_motocyles(vector<struct sub_slope> &slopes,
		vector<KLT_FeatureList> &features)
{
	const struct vehicle_mark *vm;
	int slope_width;
	int k;

	for (k = 0; k < slopes.size(); k++) {
		if (!slopes.at(k).search_plate) {
			continue;
		}
		total_slopes++;
		vm = find_vehicle_mark_in_slope(&slopes.at(k), features);
		if (NULL == vm) {
			continue;
		}
		total_slopes_gt++;
		slope_width = slopes.at(k).right - slopes.at(k).left;
		if (!vm->sema) {
			if (vm->moto) {
				slopes_motocycle_gt++;
			} else {
				slopes_car_gt++;
			}
		}
		if (slope_width > CAR_MIN_WIDTH) {
			slopes_car++;
			if (vm->moto) {
				printf("\n  Jesus, it was a motocycle...\n");
			}
		} else {
			slopes_motocycle++;
			if (!vm->moto) {
				printf("\n  God damn! Not a motocycle...\n");
			}
		}
	}
}

/**
 * Processa um array de features --- novas ou não.
 * Se o array de features estiver vazio, limpa as pirâmedes construídas.
 * Se o frame atual, cd->iframe, for o mesmo do descobrimento das features na posição k do array,
 * features.at(k)->discovered, não faz nada.
 * Se o frame atual for features.at(k)->discovered + 1 (um frame após o descobrimento), roda o SIFT
 */
void process_feature_list(vector<KLT_FeatureList> &features,
		KLT_TrackingContext tc, struct dsm_vs_data *vs,
		struct dsm_common_data *cd, const char *savelog,
		bool print, int color, uint8_t *fdebug)
{
	uint8_t *roimg;
	uint8_t *prev;
	int width;
	int height;
	int shift_x1;
	int shift_y1;
	int shift_x2;
	int shift_y2;
	int k;

	width = vs->img_width;
	height = vs->img_height;
	roimg = grab_past_buffer(cd, 0);
	prev = grab_past_buffer(cd, 1);

	if (0 == features.size()) {
		if (tc->pyramid_last != NULL) {
			_KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last);
			_KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last_gradx);
			_KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last_grady);
			tc->pyramid_last = NULL;
		}
		if (tc->pyramid_act != NULL) {
			_KLTFreePyramid((_KLT_Pyramid) tc->pyramid_act);
			_KLTFreePyramid((_KLT_Pyramid) tc->pyramid_act_gradx);
			_KLTFreePyramid((_KLT_Pyramid) tc->pyramid_act_grady);
			tc->pyramid_act = NULL;
		}
		return;
	}
	for (k = 0; k < features.size(); k++) {
		if (features.at(k)->discovered + 1 != cd->iframe) {
			continue;
		}

		/*
		printf("  DBG: rodando o SIFT: frame %d discovered %d id %d\n",
					cd->iframe, features.at(k)->discovered,
					features.at(k)->object_id);
					*/

		/**
		 * Do this only if we are in the next frame after the
		 * features selection.  Here we use the SIFT algorithm to
		 * estimate the first displacement of the vehicles.
		 */
		shift_x1 = (features.at(k)->xmin - 50 > 0 ?
				50 : features.at(k)->xmin);
		shift_y1 = (features.at(k)->ymin - 50 > 0 ?
				50 : features.at(k)->ymin);
		shift_x2 = (features.at(k)->xmin - 200 > 0 ?
				200 : features.at(k)->xmin);
		shift_y2 = (features.at(k)->ymin - 200 > 0 ?
				200 : features.at(k)->ymin);

		Image cimg1 = crop_image(prev, height, width,
				features.at(k)->xmin - shift_x1,
				features.at(k)->ymin - shift_y1,
				features.at(k)->xmax + 50,
				features.at(k)->ymax + 50);
		Keypoint key1;
		key1 = GetKeypoints(cimg1);

		Image cimg2 = crop_image(roimg, height, width,
				features.at(k)->xmin - shift_x2,
				features.at(k)->ymin - shift_y2,
				features.at(k)->xmax + 200,
				features.at(k)->ymax + 200);
		Keypoint key2;
		key2 = GetKeypoints(cimg2);

		/*
		{
			char tmpstr[256];
			sprintf(tmpstr, "/tmp/slopes/%05d_cimg1.pgm", cd->iframe);
			WritePGMFile(tmpstr, cimg1);
			sprintf(tmpstr, "/tmp/slopes/%05d_cimg2.pgm", cd->iframe);
			WritePGMFile(tmpstr, cimg2);
		}
		*/

		FindMatches(cimg1, key1, cimg2, key2, features.at(k),
				shift_x1, shift_y1, shift_x2, shift_y2,
				cd->iframe, cd->output_path);

		DisallocMatrix(cimg1->pixels, cimg1->rows,
				cimg1->cols, PERM_POOL);
		free(cimg1);
		DisallocMatrix(cimg2->pixels, cimg2->rows,
				cimg2->cols, PERM_POOL);
		free(cimg2);

		Keypoint aux;
		aux = key1;
		while (aux != NULL) {
			Keypoint next = aux->next;
			free(aux->ivec);
			free(aux);
			aux = next;
		}
		aux = key2;
		while (aux != NULL) {
			Keypoint next = aux->next;
			free(aux->ivec);
			free(aux);
			aux = next;
		}
	}

	klt_update_pyramid(features, tc, &vs->kltbufs, roimg, width, height);
	for (k = 0; k < features.size(); k++) {
		if (features.at(k)->discovered < cd->iframe) {
			/*
			printf("  DBG: rastreamento KLT: frame %d discovered %d id %d\n",
					cd->iframe, features.at(k)->discovered,
					features.at(k)->object_id);
					*/

			/* Kanade-Lucas-Tomasi Pyramidal Tracking: */
			KLTTrackFeatures(tc, width, height, features.at(k),
					cd->iframe, k);

			update_displacement(features.at(k));

			/* Speed estimation */
			compute_ipm_features(features.at(k), cd->iframe);

			outlier_removal(roimg, width, height,
					features.at(k), cd->iframe);

			speed_i speed = compute_velocity_vector(features.at(k),
					cd->iframe - 1, cd->iframe);

			features.at(k)->ipm_dx = speed.m_x -
				features.at(k)->ipm_x;
			features.at(k)->ipm_dy = speed.m_y -
				features.at(k)->ipm_y;
			features.at(k)->ipm_x = speed.m_x;
			features.at(k)->ipm_y = speed.m_y;

			features.at(k)->speed[cd->iframe % SMPSIZE] =
				metric_calc_velocity_ipm(speed, (unsigned)
						features.at(k)->vehicle.lane);
			features.at(k)->speed[cd->iframe % SMPSIZE] *= 3.6;

			if ((features.at(k)->speed[cd->iframe % SMPSIZE] < .25)
					&& (features.at(k)->discovered + 1 <
						cd->iframe) &&
					(!features.at(k)->always_track)) {
				features.at(k)->falsepositive = 1;
			}
		}
	}
	/*
	uint8_t *imgtmp;
	int size = 1.5 * width * height;
	imgtmp = (uint8_t *) malloc(size);
	memcpy(imgtmp, roimg, size);
	*/

	for (k = 0; k < features.size(); k++) {
		KLT_FeatureList kltf = features.at(k);

		if (print) {
			//vs_write_features_to_image(roimg, width, height, kltf);
			//draw_motion_vectors_nv12(roimg, width, height,
			//		cd->iframe, kltf, color);
			{
				//draw_motion_vectors_across(fdebug, width,
				//		height, cd->iframe, kltf, 4);
				//write_uchar_to_pgm(imgtmp, width, height,
				//		"/tmp/klt/%05d.pgm", cd->iframe);
			}
		}
		if (kltf->discovered + 1 < cd->iframe) {
			vs_log_vehicle_speed(kltf, cd->iframe, cd->output_path,
					savelog, &vs->should_save_frame);
			if (print && (VDSM_DLEVEL > 2)) {
				draw_speed_values(prev, width, height,
						cd->iframe, kltf,
						cd->output_path);
			}
		}
	}

	/*
	write_nv12_to_jpeg(imgtmp, width, height,
			"/tmp/klt/%05d.jpeg", cd->iframe);
	free(imgtmp);
	*/
}

void debug_print_guess_moto_stat(vector<struct sub_slope> &slopes,
		vector<KLT_FeatureList> &features)
{
#if VDSM_DLEVEL > 2
	float dmp = 0.0;
	float dcp = 0.0;

	debug_guess_motocyles(slopes, features);

	if (slopes_motocycle_gt > 0) {
		dmp = 100. * (float) slopes_motocycle /
		      (float) slopes_motocycle_gt;
	}
	if (slopes_car_gt > 0) {
		dcp = 100. * (float) slopes_car / (float) slopes_car_gt;
	}
	/*
	printf("  DBG: (%d,%d) -> "
			"MOTO(%d - %d): %.2f \tCARR(%d - %d): %.2f\n",
			total_slopes, total_slopes_gt,
			slopes_motocycle, slopes_motocycle_gt, dmp,
			slopes_car, slopes_car_gt, dcp);
			*/
#endif
}

static void process_show_speed(struct show_speed *ss,
		vector<KLT_FeatureList> f, int iframe)
{
	bool bound[NUM_SHOW_SPEED] = {false};
	int l;
	int k;

	for (k = 0; k < f.size(); k++) {
		if ((f.at(k)->discovered >= iframe) ||
				(!f.at(k)->vehicle.radar) ||
				(f.at(k)->vehicle.lane < 1) ||
				(f.at(k)->vehicle.lane > 3)) {
			continue;
		}
		l = f.at(k)->vehicle.lane - 1;
		if ((f.at(k)->ymin < 320) && (!f.at(k)->print_speed)) {
			f.at(k)->print_speed = true;
			ss[l].wait_meas = true;
			ss[l].cnt_show = 120;
			ss[l].speed = f.at(k)->vehicle.speed;
			ss[l].show_m = false;
			ss[l].sum_meas = 0.0;
			ss[l].cnt_meas = 0;
		}
		if (f.at(k)->print_speed && (f.at(k)->ymin < 320)) {
			bound[l] = true;
			ss[l].sum_meas += f.at(k)->speed[iframe % SMPSIZE];
			ss[l].cnt_meas++;
		}
	}
	for (k = 0; k < NUM_SHOW_SPEED; k++) {
		if (!ss[k].wait_meas) {
			continue;
		}
		if (!bound[k] && (ss[k].cnt_meas > 0)) {
			ss[k].measured =
				ss[k].sum_meas / (float) ss[k].cnt_meas;
			ss[k].show_m = true;
		}
	}
}

static void print_vehicles_speed(uint8_t *img, int w, int h,
		struct show_speed *ss)
{
	int cb;
	int cr;
	int k;

	for (k = 0; k < NUM_SHOW_SPEED; k++) {
		if (ss[k].cnt_show <= 0) {
			continue;
		}
		ss[k].cnt_show--;
		draw_text(img, w, ss[k].x, ss[k].y, 4, 0, 220,
				"   Speed:  ");
		draw_text(img, w, ss[k].x, ss[k].y + 38, 4, 0, 220,
				" %.1f km/h ", ss[k].speed);
		if (ss[k].show_m) {
			//if (ss[k].measured - ss[k].speed > 5) {
			//	ss[k].measured = ss[k].speed - 0.6;
			//}
			draw_text(img, w, ss[k].x, ss[k].y + 86, 4, 0, 220,
					" Measured: ");
			draw_text(img, w, ss[k].x, ss[k].y + 124, 4, 0, 220,
					" %.1f km/h ", ss[k].measured);
			if ((ss[k].measured - ss[k].speed) > 2.0) {
				cb = 82;
				cr = 240;
			} else if ((ss[k].measured - ss[k].speed) < - 3.0) {
				cb = 240;
				cr = 24;
			} else {
				cb = 26;
				cr = 26;
			}
			paint_chrom_rect(img, w, h, ss[k].x, ss[k].y + 86,
					ss[k].x + 352, ss[k].y + 160, cb, cr);
		}
	}
}

/**
 * Main callback function.
 * @param buff Pointer to current image buffer.
 * @param vs Pointer to the structure returned by vs_preload function.
 */
void vs_chain_callback(struct dsm_common_data *cd, void *_data)
{
	static vector<struct sub_slope> slopes;
	vector<plate> plates_gt;
	struct dsm_vs_data *vs;
	struct sub_data *sub;
	uint8_t *roimg;
	uint8_t *prev;
	uint8_t *chrom;
	uint8_t *fdebug;
	int width;
	int height;
	bool ret = 0;
	int k;
	static int shift_iframe = 0;
	int hold_slope = 0;

	/* Get structs. */
	vs = (struct dsm_vs_data *) _data;
	width = vs->img_width;
	height = vs->img_height;
	sub = vs->background_sub;

	/* Get buffers. */
	roimg = grab_past_buffer(cd, 0);
	prev = grab_past_buffer(cd, 1);
	chrom = roimg + (width * height);
	decolorize(roimg, width, height);

	print_timestamp("frame " cterm(TB_WHT, "%d\n"), cd->iframe);
	/*
	if (cd->iframe < 430) {
		return;
	}
	*/

	fdebug = (unsigned char *) malloc(width * height * 1.5);
	memcpy(fdebug, roimg, width * height * 1.5);

	plates_gt = read_plates_from_xml(vs->vt, cd->iframe);
	if (vs->vt && (cd->iframe > vs->vt->size)) {
		print_warn("Reached end of ground truth file in frame %d\n",
				cd->iframe);
		vs_end_of_stream_handler(_data);
	}

	/* Start from here doing the background subtraction. */
	ret = background_subtraction(sub, roimg, cd->iframe, cd->timestamp,
			fdebug);
	if (ret && use_detection) {
		/* If there is some significant changes in background. */
		slopes = compute_region_of_interest(sub, roimg, slopes);

		/**
		 * Set the flag area[lane].tracking_plate if there is a
		 * license plate being tracked in the corresponding lane.
		 */
		print_tab(DBG_CALL, "Check lane already tracking\n");
		check_lane_already_tracking(&vs->lanes, vs->features,
				MINY_TRACKING);
		print_tab(DBG_EXIT, "Check lane already tracking\n");

		/* Set or reset the flag search_plate for each slope. */
		print_tab(DBG_CALL, "Update search plate\n");
		update_slope_search_plate(slopes, &vs->lanes, height -
				WINDOW_BOTTON_MARGIN - SUBSAMPLE_STEP_Y);
		print_tab(DBG_EXIT, "Update search plate\n");
		check_slope_already_tracking(vs->features, slopes);
	} else {
		if (slopes.size() > 0) {
			slopes.clear();
		}
		for (k = 0; k < sub->hprof.size; k++) {
			sub->hprof.profile[k] = 0;
		}
	}
	sub_paint_debug(sub, fdebug, slopes);
	/* Remove slopes with flag search_plate = false. */
	for (k = 0; k < slopes.size();) {
		if (!slopes.at(k).search_plate) {
			slopes.erase(slopes.begin() + k);
		} else {
			k++;
		}
	}

	if (slopes.size() > 0) {
		vector<plate> plates;
#if 0
		print_tab(DBG_CALL, "SnooperText\n");
		plates = snoopertext_detection_working(
				roimg,
				&vs->snooper,
				height,
				width,
				cd->iframe,
				SUBSAMPLE_STEP_X,
				slopes,
				WINDOW_SHIFT_X_MAX,
				WINDOW_SHIFT_Y_MAX);
		print_tab(DBG_EXIT, "SnooperText\n");
#endif
#if 1
		print_tab(DBG_CALL, "Detection Minetto\n");
		char imgname[32];
		sprintf(imgname, "%05d", cd->iframe);
		plates = detection_minetto(
				roimg,
				&vs->snooper,
				width,
				height,
				(char *) "/tmp/slopes",
				imgname,
				slopes,
				WINDOW_SHIFT_X_MAX,
				WINDOW_SHIFT_Y_MAX,
				cd->iframe,
                                sub->umask);
		print_tab(DBG_EXIT, "Detection Minetto\n");
#endif
#if 0
		print_tab(DBG_CALL, "Detection Zheng\n");
		char imgname[32];
		sprintf(imgname, "%05d", cd->iframe);
		plates = detection_zheng(
				roimg,
				&vs->snooper,
				width,
				height,
				(char *) "/tmp/slopes",
				imgname,
				slopes,
				WINDOW_SHIFT_X_MAX,
				WINDOW_SHIFT_Y_MAX,
				cd->iframe,
                                sub->umask);
		print_tab(DBG_EXIT, "Detection Zheng\n");
#endif
#if 0
		print_tab(DBG_CALL, "Detection SWT\n");
		plates = detection_swt(
				roimg,
				&vs->snooper,
				width,
				height,
				(char *) "/tmp/slopes",
				slopes,
				WINDOW_SHIFT_X_MAX,
				WINDOW_SHIFT_Y_MAX,
				cd->iframe);
		print_tab(DBG_EXIT, "Detection SWT\n");
#endif

#if VDSM_DLEVEL > 1
		dsm_paint_plates(roimg, width, height, plates, 255, 0);
#endif
		print_tab(DBG_CALL, "Select Plates\n");
		dsm_select_plates(plates, &vs->lanes, cd->iframe);
		print_tab(DBG_EXIT, "Select Plates\n");
		print_tab(DBG_CALL, "Search ground truth from XML\n");
		search_groundtruth_from_xml(plates, vs->vt, cd->iframe);
		print_tab(DBG_EXIT, "Search ground truth from XML\n");
#if VDSM_DLEVEL > 1
		sub_paint_debug(sub, roimg, slopes);
#endif
#if VDSM_DLEVEL > 1
		dsm_paint_plates(fdebug, width, height, plates, 72, 220);
#endif
#if 0
		save_slopes_to_pgm(roimg, width, height, "/tmp",
				slopes, cd->iframe);
#endif
		print_tab(DBG_CALL, "Compute features from plates\n");
		compute_features_from_plates(plates, vs->features, vs, cd,
				false);
		print_tab(DBG_EXIT, "Compute features from plates\n");

		/*
		if (plates.size() > 0) {
			hold_slope = 15;
		}
		*/
		plates.clear();
	}
	//compute_features_from_plates(plates_gt, vs->features_gt, vs, cd, true);
	//dsm_paint_plates(roimg, width, height, plates_gt, 0, 255);

	debug_print_guess_moto_stat(slopes, vs->features_gt);

	vs_print_vehicle_plates(roimg, vs, cd->iframe, 5);
	vs_unnecessary_remove_features(vs->features);
	vs_unnecessary_remove_features(vs->features_gt);

	print_tab(DBG_CALL, "Process features list\n");
	process_feature_list(vs->features, vs->tc, vs, cd,
			"detection", true, 255, roimg);
	process_feature_list(vs->features_gt, vs->tc_gt, vs, cd,
			"groundtruth", true, 0, fdebug);
	print_tab(DBG_EXIT, "Process features list\n");

	cd->timestamp += VIDEO_FRAME_PERIOD;

	/*
	for (k = 0; k < slopes.size();) {
		if (slopes.at(k).search_plate) {
			hold_slope = 15;
			break;
		}
	}
	*/
	process_show_speed(sspeed, vs->features, cd->iframe);
	print_vehicles_speed(roimg, width, height, sspeed);
	/*
	write_nv12_to_jpeg(fdebug, width, height, "/tmp/fdebug/%05d.jpeg",
			cd->iframe + shift_iframe);
	while (hold_slope--) {
		shift_iframe++;
		write_nv12_to_jpeg(fdebug, width, height,
				"/tmp/fdebug/%05d.jpeg",
				cd->iframe + shift_iframe);
	}
	*/
	free(fdebug);

	/* Paint the chrominance from image, if enabled. */
#if VDSM_DLEVEL > 2
	vs_print_lanes_to_nv12(chrom, width, height, vs->lanes.lmask, 0);
	vs_print_lanes_to_nv12(chrom, width, height, vs->lanes.mmask, 32);
	if (vs->should_save_frame) {
		print_dbg("saving frame %d\n", cd->iframe);
		write_nv12_to_jpeg(roimg, width, height, "%s/snapshot/%05d.jpg",
				cd->output_path, cd->iframe);
	}
#endif
}

}
