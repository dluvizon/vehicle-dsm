/**
 * file: backgroundsub.cpp
 */

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "backgroundsub.h"
#include "straux.h"
#include "fsaux.h"
#include "font.h"
#include "image.h"
#include "utils.h"
#include "config.h"

#include "sift/key.h"

//#include "opencv_canny.h"

#define _DLEVEL BACKGROUNDSUB_DLEVEL
#include "debug.h"

static void alloc_profile(struct profile *prof, int size)
{
	prof->profile = (int *) malloc(size * sizeof(int));
	prof->min = (int *) malloc(size * sizeof(int));
	prof->max = (int *) malloc(size * sizeof(int));
	prof->smooth = (int *) malloc(size * sizeof(int));
	prof->rising = (int *) malloc(size * sizeof(int));
	prof->falling = (int *) malloc(size * sizeof(int));
	prof->accum = (int *) calloc(sizeof(int), size * sizeof(int));
	prof->size = size;
}

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

/**
 * Configure the window to compute background subtraction.
 */
static void sub_set_window(struct sub_data *sdata, struct rectangle_m window,
		int xstep, int ystep)
{
	int dx, dy;
	int size;

	sdata->win = window;
	sdata->xstep = xstep;
	sdata->ystep = ystep;
	dx = window.xmax - window.xmin;
	dy = window.ymax - window.ymin;
	assert(dx > 0 && dy > 0);
	sdata->area =  dx * dy / (xstep * ystep);

	/* Alloc memory for the struct horizontal profile. */
	size = (sdata->width / xstep) + 1;
	alloc_profile(&sdata->hprof, size);

	/* Alloc memory for the struct vertical profile. */
	size = (sdata->height / ystep) + 1;
	alloc_profile(&sdata->vprof, size);
}

struct sub_data *background_sub_new(int width, int height,
		struct rectangle_m window, int xstep, int ystep)
{
	struct sub_data *sdata;
	int i;

	sdata = (struct sub_data *) malloc(sizeof(struct sub_data));
	if (!sdata) {
		print_err("could not alloc memory for sub_data");
		goto end;
	}
	memset(sdata, 0, sizeof(struct sub_data));
	sdata->mhi = (double *) malloc(width * height * sizeof(double));
	sdata->umask = (uint8_t *) malloc(width * height * sizeof(uint8_t));
	sdata->buffers = (double **) malloc(SUB_NFRAMES * sizeof(double *));
	for (i = 0; i < SUB_NFRAMES; i++) {
		sdata->buffers[i] = (double *)
			malloc(width * height * sizeof(double));
		fill_double_image(sdata->buffers[i], height, width, 0.0);
	}
	sdata->width = width;
	sdata->height = height;
	sdata->threshold = SUB_THRESHOLD;

	sub_set_window(sdata, window, xstep, ystep);
	sdata->slope_id_cnt = 0;

	/* Configure max and min for each lane. */
	sdata->xlane_min[0] = 10 / xstep;
	sdata->xlane_max[0] = 565 / xstep;
	sdata->xlane_min[1] = 505 / xstep;
	sdata->xlane_max[1] = 1250 / xstep;
	sdata->xlane_min[2] = 1200 / xstep;
	sdata->xlane_max[2] = 1880 / xstep;

end:
	return sdata;
}

void sub_create_logdir(struct sub_data *sdata, const char *path)
{
	size_t len;
	int i;

	assert(sdata && path);

	/* first check if the directory exist, if not, create it */
	if (!fs_dir_exist(path)) {
		/* this directory does not exist */
		print_warn("attempt to create directory %s\n", path);
		if (fs_mkdir_recursive(path)) {
			print_err("fail to create dir!");
			return;
		}
	}

	/* now copy the path string to a local variable */
	len = strlen(path);
	sdata->logdir_path = (char *) malloc(len + 1);
	if (!sdata->logdir_path) {
		print_err("could not alloc memory for logdir_path");
		return;
	}
	strcpy(sdata->logdir_path, path);
	for (i = 0; i < GTRUTH_MAX_LANES; i++) {
		char tmp[256];

		snprintf(tmp, sizeof(tmp), "%s/stat_f%d.log", path, i + 1);
		sdata->stat[i].logfile = fopen(tmp, "w");
	}
}

void sub_write_stat_log(struct sub_stat_data *stat)
{
	float err = 0.0;
	int ret;

	if (!stat->logfile) {
		print_err("log file NULL");
		return;
	}
	ret = fseek(stat->logfile, 0L, SEEK_SET);
	if (-1 == ret) {
		print_err("error resetting stat file");
		return;
	}
	if (stat->cnt_vehicles > 0) {
		err = (float) stat->cnt_slopes / (float) stat->cnt_vehicles -1.;
	}
	fprintf(stat->logfile, "vehicles   %ld\n", stat->cnt_vehicles);
	fprintf(stat->logfile, "background %ld\n", stat->cnt_slopes);
	fprintf(stat->logfile, "err        %.2f %%\n", err * 100.0);
	if (stat->cnt_flush++ >= 10) {
		stat->cnt_flush = 0;
		fflush(stat->logfile);
	}
}

static void sub_paint_pixel(uint8_t *img, int w, int h, int x, int y, bool bold)
{
	uint8_t *chrom;
	uint8_t cb = 255;
	uint8_t cr = 0;
	uint8_t ym = 128;
	int pw = 8;
	int ph = 6;
	int iw;
	int ih;

	if (x < pw) {
		x = pw;
	}
	if (y < ph) {
		y = ph;
	}
	if (x > w - pw - 1) {
		x = w - pw - 1;
	}
	if (y > h - ph - 1) {
		y = h - ph - 1;
	}

	chrom = img + w * h;
	for (iw = - pw / 2; iw <= pw / 2; iw++) {
		for (ih = - ph / 2; ih <= ph / 2; ih++) {
			print_pixel_lum(img, w, h, x + iw, y + ih, ym);
			print_pixel_nv12(chrom, w, h, x + iw, y + ih, cb, cr);
		}
	}
	/*
	print_pixel_nv12(chrom, w, h, x, y, 255, 0);
	if (bold) {
		print_pixel_nv12(chrom, w, h, x - 2, y,     255, 0);
		print_pixel_nv12(chrom, w, h, x + 2, y,     255, 0);
		print_pixel_nv12(chrom, w, h, x,     y - 2, 255, 0);
		print_pixel_nv12(chrom, w, h, x - 2, y - 2, 255, 0);
		print_pixel_nv12(chrom, w, h, x + 2, y - 2, 255, 0);
		print_pixel_nv12(chrom, w, h, x,     y + 2, 255, 0);
		print_pixel_nv12(chrom, w, h, x - 2, y + 2, 255, 0);
		print_pixel_nv12(chrom, w, h, x + 2, y + 2, 255, 0);
	}
	*/
}

static void plot_vertical_profile(uint8_t *img, int width, int height,
		int ystep, int offset, struct profile *prof, bool bold)
{
	int x;
	int y;
	int i;

	for (i = 2; i < prof->size - 2; i++) {
		y = i * ystep;
		x = offset + abs(prof->profile[i] / SUBSAMPLE_STEP_X);
		if (x < 2) {
			x = 2;
		}
		if (x > width - 3) {
			x = width - 3;
		}
		if (y < 2) {
			y = 2;
		}
		if (y > height - 3) {
			y = height - 3;
		}
		sub_paint_pixel(img, width, height, x, y, bold);
	}
}

/**
 * Plot the profile stored in prof->profile at the yplot height.
 * TODO: params
 */
static void sub_plot_profile(uint8_t *img, int width, int height, int xstep,
		struct profile *prof, bool bold)
{
	int x, y;
	int i;

	/* plot the horizontal profile in the middle of the image */
	for (i = 2, x = 0; i < prof->size - 2; i++, x += xstep) {
		y = 880 - abs(prof->profile[i] / SUBSAMPLE_STEP_Y);
		sub_paint_pixel(img, width, height, x, y, bold);
	}
}


/**
 * If the umask is significant, update the window of events and return 1.
 * Otherwise, return zero.
 */
static int sub_detecting_changes(struct sub_data *sdata, int i, int x, int y)
{
	if (sdata->umask[i] > 127) {
		if (x < sdata->wevent.xmin) {
			sdata->wevent.xmin = x;
		}
		if (x > sdata->wevent.xmax) {
			sdata->wevent.xmax = x;
		}
		if (y < sdata->wevent.ymin) {
			sdata->wevent.ymin = y;
		}
		if (y > sdata->wevent.ymax) {
			sdata->wevent.ymax = y;
		}
		return 1;
	}
	return 0;
}

bool background_subtraction(struct sub_data *sdata, uint8_t *frame, int iframe,
		double timestamp, uint8_t *fdebug)
{
	struct rectangle_m *win;
	double *silhouette;
	int motion = 0;
	int ifl;
	int ifx;
	int x;
	int y;

	/*
	 * Index of the (last - (SUB_NFRAMES-1))th frame.
	 * This will store the current frame.
	 */
	ifl = sdata->iframe_last;
	/* Index of the SUB_NFRAMES th frame, pointed by silhouette. */
	ifx = (ifl + 1) % SUB_NFRAMES;
	silhouette = sdata->buffers[ifx];

	/* Reset the window of events. */
	sdata->wevent.xmin = sdata->width;
	sdata->wevent.xmax = 0;
	sdata->wevent.ymin = sdata->height;
	sdata->wevent.ymax = 0;

	/* Run background subtraction only inside the window area. */
	win = &sdata->win;

	/* Iterate inside the window area, from x to y. */
	for (y = win->ymin; y < win->ymax; y += sdata->ystep) {
		for (x = win->xmin; x < win->xmax; x += sdata->xstep) {
			int i = y * sdata->width + x;

			/* Copying the original image to the double buffer. */
			sdata->buffers[ifl][i] = (double) (frame[i]);

			/*
			 * Adaptive background subtraction (absolute frame
			 * difference pixel-by-pixel).
			 */
			silhouette[i] = abs(sdata->buffers[ifl][i] -
					sdata->buffers[ifx][i]);

			/*
			 * Updating the {mhi} image, which contains the video
			 * motion history (the background model).
			 */
			if (silhouette[i] >= sdata->threshold) {
				/* New motions detected. */
				sdata->mhi[i] = timestamp;
			}
			if (sdata->mhi[i] > (timestamp - SUB_MHI_DURATION)) {
				/* Changed before the MHI duration timeout. */
				sdata->umask[i] = 255;
			} else {
				/* Changes too old. */
				sdata->umask[i] = 0;
			}
			motion += sub_detecting_changes(sdata, i, x, y);
		}
	}

	/* Changing the indice of the last frame. */
	sdata->iframe_last = ifx;

	/* Write log images if enabled in condig.h. */
#if BACKGROUNDSUB_DLEVEL > 3
	write_double_to_pgm(silhouette, sdata->width, sdata->height,
			"/tmp/slopes/sil_%05d.pgm", iframe);
#endif

#if 0
	{
		char fname[256];
		unsigned char *imgchr;

		imgchr = fdebug + (sdata->width * sdata->height);
		for (y = win->ymin; y < win->ymax; y += 1) {
			for (x = win->xmin; x < win->xmax; x += 1) {
				fdebug[y * sdata->width + x] = 0;
				print_pixel_nv12(imgchr, sdata->width,
						sdata->height, x, y, 128, 128);
			}
		}
		for (y = win->ymin; y < win->ymax; y += sdata->ystep) {
			for (x = win->xmin; x < win->xmax; x += sdata->xstep) {
				int i = y * sdata->width + x;
				if (sdata->umask[i]) {
					fdebug[i] = 255;
				}
			}
		}
	}
#endif

	/* Write debug image to output video stream if enabled in config.h. */
#if BACKGROUNDSUB_DLEVEL > 2
	memcpy(frame, sdata->umask, sdata->width * sdata->height);
	draw_rectanglem_to_uchar(frame, sdata->width, &sdata->wevent);
#endif

	if (((double) motion / sdata->area) > SUB_PERCENTAGE) {
		return true;
	}
	return false;
}

static int check_matching_id(struct sub_slope &slope,
		vector<struct sub_slope> &prev_slopes)
{
	int i;

	for (i = 0; i < prev_slopes.size(); i++) {
		if ((slope.right < prev_slopes.at(i).left) ||
				(slope.left > prev_slopes.at(i).right)) {
			continue;
		}
		slope.id = prev_slopes.at(i).id;
                slope.ntries = prev_slopes.at(i).ntries + 1;
		prev_slopes.erase(prev_slopes.begin() + i);
		break;
	}
}

static void fit_slope_inside_wevent(struct sub_slope *slope,
		struct rectangle_m *wevent)
{
	if (slope->up < wevent->ymin) {
		slope->up = wevent->ymin;
	}
	if (slope->down > wevent->ymax) {
		slope->down = wevent->ymax;
	}
	if (slope->left < wevent->xmin) {
		slope->left = wevent->xmin;
	}
	if (slope->right > wevent->xmax) {
		slope->right = wevent->xmax;
	}
}

/**
 * TODO: understand in a better way this filter.
 */
/**
 * Smooth the profile curve.
 * @param prof Pointer to a horizontal profile.
 * @param iter Number of iterations to run the algorithm.
 * @param wsize Window size.
 */
static void window_smooth(struct profile *prof, int iter, int wsize)
{
	double sum;
	int i, j;

	int_array_set(prof->smooth, 0, prof->size);
	while (iter--) {
		for (j = wsize / 2; j < prof->size - wsize / 2; j++) {
			sum = 0.0;
			for (i = -wsize / 2; i <= wsize / 2; i++) {
				sum += prof->profile[j + i];
			}
			prof->smooth[j] = (int)((double) sum / (double) wsize);
		}
		for (j = 0; j < prof->size; j++) {
			if ((j < wsize / 2) || (j > prof->size - wsize / 2)) {
				prof->profile[j] = 0;
			} else {
				prof->profile[j] = prof->smooth[j];
			}
		}
	}
}

/**
 * This function returns a vector with the values with the time since
 * the last rise/fall, according to the r vector.
 */
int *funcA(int *r, int sizer, int size, int *dr)
{
	int *rc;
	int i;

	dr[1] = r[1];
	for (i = 2; i < sizer; i++) {
		dr[i] = r[i] - r[i - 1];
	}

	/* repmat (1, nx, 1) where nx = size */
	rc = (int *) malloc(size * sizeof(int));
	for (i = 1; i < size; i++) {
		rc[i] = 1;
	}
	for (i = 2; i <= sizer; i++) {
		rc[r[i - 1] + 1] = 1 - dr[i - 1];
	}
	rc[0] = 0;
	rc[1] = 0;

	int *rs = (int *) malloc(size * sizeof(int));
	for (i = 1; i < size; i++) {
		rs[i] = 0;
	}

	for (i = 2; i < size; i++) {
		rs[i] = rc[i] + rs[i - 1];
	}
	free(rc);

	return rs;
}

int *funcB (int *r, int *dr, int sizer, int size) {
   int i;
   int *rp = (int *)malloc(size * sizeof(int));
   for (i = 1; i < size; i++) { rp[i] = -1; }
   rp[1] = dr[1] - 1;
   for (i = 1; i < sizer-1; i++) { rp[r[i]+1] = dr[i+1] - 1; }
   rp[r[i]+1] = (size-1) - r[sizer-1] - 1;
   int *rq = (int *)malloc(size * sizeof(int));
   rq[1] = rp[1];
   for (i = 2; i < size; i++) { rq[i] = rp[i] + rq[i-1]; }
   free(rp);
   return rq;
}

static void reset_sub_slopes(struct sub_slope *s)
{
	s->mag = 0;
	s->up = 0;
	s->down = 0;
	s->left = 0;
	s->center = 0;
	s->right = 0;
	s->id = 0;
	s->search_plate = false;
}

vector<struct sub_slope> findpeaks(struct profile *prof, int size)
{
	vector<struct sub_slope> slopes;
	int sizer = 1;
	int sizef = 1;
	int *rising;
	int *falling;
	int delta;
	int *v;
	int i;

	v = prof->profile;
	rising = prof->rising;
	falling = prof->falling;

	/**
	 * Find the rising and falling borders in profile, which is
	 * identified by a rising or a falling changing in value
	 * greater than PROF_DELTA_PERC.
	 */
	for (i = 1; i < size - 1; i++) {
		delta = v[i + 1] - v[i];
		if (v[i] != 0) {
			delta = (100 * delta) / v[i];
		} else if (delta > 0) {
			delta = 100;
		} else if (delta < 0){
			delta = -100;
		}
		if (delta > PROF_DELTA_PERC) {
			rising[sizer++] = i;
		} else if (delta < -PROF_DELTA_PERC) {
			falling[sizef++] = i;
		}
	}

	/* If we found at least one rising and one falling border. */
	if ((sizer > 1) && (sizef > 1)) {
		int *dr = (int *) malloc(sizer * sizeof(int));
		int *rs = funcA(rising, sizer, size, dr);
		int *rq = funcB(rising, dr, sizer, size);

		int *df = (int *) malloc(sizef * sizeof(int));
		int *fs = funcA(falling, sizef, size, df);
		int *fq = funcB(falling, df, sizef, size);

		int *ds = (int *) malloc(size * sizeof(int));
		int *dq = (int *) malloc(size * sizeof(int));
		int *dt = (int *) malloc(size * sizeof(int));

		for (i = 0; i < size; i++) {
			ds[i] = (rs[i] < fs[i] ? 1 : 0);
			dq[i] = (fq[i] < rq[i] ? 1 : 0);
			dt[i] = ((int) (floor((fq[i] - rs[i]) / 2.0)) == 0 ?
					1 : 0);
		}

		free(fs);
		free(fq);
		free(dr);
		free(rs);
		free(rq);
		free(df);

		int pleft = 0;
		int pcenter = 0;
		int pright = 0;

		for (i = 1; i < size - 1; i++) {
			/* Ascending part: */
			if (!ds[i - 1] && ds[i]) {
				struct sub_slope s;
				reset_sub_slopes(&s);
				s.left = i - 1;
				slopes.push_back(s);
				pleft++;
			}
			/* Peak center: */
			if (ds[i] && dq[i] && dt[i]) {
				slopes.at(pcenter++).center = i - 1;
			}
			/* Descending part: */
			if ((dq[i] && !dq[i + 1]) && (pleft > pright)) {
				slopes.at(pright++).right = i - 1;
			}
		}
		free(ds);
		free(dq);
		free(dt);
	}

	return slopes;
}

static int get_nearest_ilane(int center, int *vmid)
{
	int min_dist = INT_MAX;
	int ilane = 0;
	int k;

	for (k = 0; k < 3; k++) {
		if (min_dist > abs(center - vmid[k])) {
			min_dist = abs(center - vmid[k]);
			ilane = k;
		}
	}

	return ilane;
}

/**
 * Split the slopes according to the lane's boundary.
 */
static void split_slopes(vector<struct sub_slope> &slopes,
		int *xlane_max, int *xlane_min)
{
	struct sub_slope s;
	int mid[3];
	int ilane;
	int i, k;

	for (k = 0; k < 3; k++) {
		mid[k] = (xlane_max[k] + xlane_min[k]) / 2;
	}

	for (i = 0; i < slopes.size();) {
		ilane = get_nearest_ilane(slopes.at(i).center, mid);
		if (slopes.at(i).left < xlane_min[ilane]) {
			/* Split this in the left side. */
			if (ilane > 0) {
				/* Create new in the left.*/
				s = slopes.at(i);
				s.right = xlane_max[ilane - 1];
				if ((s.right - s.left) >= (2 * (slopes.at(i).right - xlane_min[ilane]) / 3)) {
					s.center = (s.left + s.right) / 2;
					if (get_nearest_ilane(s.center, mid) !=
							ilane) {
						slopes.push_back(s);
					}
				}
			}
			slopes.at(i).left = xlane_min[ilane];
			slopes.at(i).center =
				(slopes.at(i).left + slopes.at(i).right) / 2;
		}
		if (slopes.at(i).right > xlane_max[ilane]) {
			/* Split this in the right side. */
			if (ilane < 2) {
				/* Create new in the right.*/
				s = slopes.at(i);
				s.left = xlane_min[ilane + 1];
				if ((s.right - s.left) > (2 * (xlane_max[ilane] - slopes.at(i).left) / 3)) {
					s.center = (s.left + s.right) / 2;
					if (get_nearest_ilane(s.center, mid) !=
							ilane) {
						slopes.push_back(s);
					}
				}
			}
			slopes.at(i).right = xlane_max[ilane];
			slopes.at(i).center =
				(slopes.at(i).left + slopes.at(i).right) / 2;
		}
		i++;
	}
}

static void compute_horizontal_profile(struct sub_data *sdata)
{
	struct profile *prof;
	struct rectangle_m *win;
	uint8_t *grid;
	int width;
	int xstep;
	int ystep;
	int x, y;
	int p;

	prof = &sdata->hprof;
	win = &sdata->win;
	grid = sdata->umask;
	width = sdata->width;
	xstep = sdata->xstep;
	ystep = sdata->ystep;

	int_array_set(prof->min, INT_MAX, prof->size);
	int_array_set(prof->max, 0, prof->size);
	for (x = win->xmin; x < win->xmax; x += xstep) {
		p = x / xstep + 1;
		prof->accum[p] = 0;
		for (y = win->ymin; y < win->ymax; y += ystep) {
			if (grid[y * width + x] > 0) {
				prof->accum[p] += SUBSAMPLE_AREA;
				if (y < prof->min[p]) {
					prof->min[p] = y;
				}
				if (y > prof->max[p]) {
					prof->max[p] = y;
				}
			}
		}
	}
	memcpy(prof->profile, prof->accum, prof->size * sizeof(int));
}

static void compute_vertical_profile(struct sub_data *sdata,
		int xmin, int xmax, int ymin, int ymax)
{
	struct profile *prof;
	uint8_t *grid;
	int width;
	int xstep;
	int ystep;
	int x, y;
	int p;

	prof = &sdata->vprof;
	grid = sdata->umask;
	width = sdata->width;
	xstep = sdata->xstep;
	ystep = sdata->ystep;

	int_array_set(prof->min, INT_MAX, prof->size);
	int_array_set(prof->max, 0, prof->size);
	int_array_set(prof->accum, 0, prof->size);

	for (y = ymin; y < ymax; y += ystep) {
		p = y / ystep + 1;
		for (x = xmin; x < xmax; x += xstep) {
			if (grid[y * width + x] > 0) {
				prof->accum[p] += SUBSAMPLE_AREA;
				if (x < prof->min[p]) {
					prof->min[p] = x;
				}
				if (x > prof->max[p]) {
					prof->max[p] = x;
				}
			}
		}
	}
	memcpy(prof->profile, prof->accum, prof->size * sizeof(int));
}

static void edge_vertical_profile(struct profile *prof, uint8_t *edge,
		int width, int height, int step, int xmin, int ymin)
{
	int x, y;
	int p;

	int_array_set(prof->min, INT_MAX, prof->size);
	int_array_set(prof->max, 0, prof->size);
	int_array_set(prof->accum, 0, prof->size);

	for (y = 0; y < height; y += 1) {
		p = (y + ymin) / step;
		for (x = 0; x < width; x += 1) {
			if (edge[x + y * width] > 0) {
				prof->accum[p] += 10;
				if ((x + xmin) < prof->min[p]) {
					prof->min[p] = (x + xmin);
				}
				if ((x + xmin) > prof->max[p]) {
					prof->max[p] = (x + xmin);
				}
			}
		}
	}
	memcpy(prof->profile, prof->accum, prof->size * sizeof(int));
}

static void find_up_down(struct sub_slope *s, struct profile *prof, int step)
{
	int th;
	int min;
	int max;
	int vmax;
	int imax;
	int i;
	bool found;

	min = s->up;
	max = s->down;
	vmax = 0;

	th = 0;
	for (i = min; i <= max; i += step) {
		if (prof->profile[i / step] > vmax) {
			vmax = prof->profile[i / step];
			imax = i;
		}
		th += prof->profile[i / step];
	}
	th /= (max - min) / step;
	th = (int) (0.30 * (float) th);
	if (th < 8 * step) {
		th = 8 * step;
	}

	for (i = imax; i >= min; i -= step) {
		if (prof->profile[i / step] < th) {
			s->up = i;
			break;
		}
	}
	found = false;
	for (i = imax; i <= max; i += step) {
		if (!found) {
			if (prof->profile[i / step] < th) {
				s->down = i;
				found = true;
			}
		} else {
			if (prof->profile[i / step] > th) {
				found = false;
			}
		}
	}
}

static void sobel_win(uint8_t *edge, int we, uint8_t *src, int ws,
		int xstep, int ystep, int thr)
{
	int x, y;
	int pix1;
	int pix2;
	int pe;
	int ps;

	for (y = 0; y < ystep; y++) {
		for (x = 0; x < xstep; x++) {
			ps = x + y * ws;
			pe = x + y * we;

			pix1 = -(int) src[ps - 1 - ws];
			pix1 -= 2 * (int) src[ps - 1 + 0];
			pix1 -= (int) src[ps - 1 + ws];
			pix1 += (int) src[ps + 1 - ws];
			pix1 += 2 * (int) src[ps + 1 + 0];
			pix1 += (int) src[ps + 1 + ws];
			pix1 = abs(pix1);

			pix2 = -(int) src[ps - 1 - ws];
			pix2 -= 2 * (int) src[ps + 0 - ws];
			pix2 -= (int) src[ps + 1 - ws];
			pix2 += (int) src[ps - 1 + ws];
			pix2 += 2 * (int) src[ps + 0 + ws];
			pix2 += (int) src[ps + 1 + ws];
			pix2 = abs(pix2) + pix1;
			if (pix2 > thr) {
				edge[pe] = 255;
			} else {
				edge[pe] = 0;
			}
		}
	}
}

static void save_slope(uint8_t *frame, int w, int h,
		int x1, int y1, int x2, int y2, int cnt)
{
	uint8_t *aux;
	uint8_t *ptr;
	uint8_t *src;
	int tw;
	int th;
	int y;

	tw = x2 - x1;
	th = y2 - y1;
	aux = new uint8_t[tw * th];
	ptr = aux;
	src = frame + x1 + y1 * w;
	for (y = y1; y < y2; y++) {
		memcpy(ptr, src, tw);
		ptr += tw;
		src += w;
	}
	write_uchar_to_pgm(aux, tw, th, "/tmp/%05d-frame.pgm", cnt);

	delete aux;
}

static void copy_sift_image_to_uchar(uint8_t *dest, Image img, int w, int h)
{
	float max = FLT_MIN;
	float min = FLT_MAX;
	float slope;
	int x, y;
	int p;

	assert((w == img->cols) && (h == img->rows));

	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			if (img->pixels[y][x] > max) {
				max = img->pixels[y][x];
			}
			if (img->pixels[y][x] < min) {
				min = img->pixels[y][x];
			}
		}
	}
	slope = 255.0 / (max - min);

	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			p = x + y * w;
			dest[p] = (uint8_t) (slope * (img->pixels[y][x] - min));
		}
	}
}

static int test_slope_up_and_down(struct sub_slope *s, uint8_t *frame,
		struct sub_data *sdata)
{
	Image image;
	static int cnt = 0;
	uint8_t *smooth;
	uint8_t *edge;
	int pdest;
	int psmooth;
	int psrc;
	int ydest;
	int xdest;
	int ysrc;
	int xsrc;
	int w;
	int h;

	image = crop_image(frame, sdata->height, sdata->width,
			s->left - 1, s->up - 1, s->right, s->down);
	GaussianBlur(image, 2);

	w = s->right - s->left;
	h = s->down - s->up;

	smooth = new uint8_t[(w + 2) * (h + 2)];
	edge = new uint8_t[w * h];
	memset(edge, 0, w * h * sizeof(uint8_t));
	memset(smooth, 0, (w + 2) * (h + 2) * sizeof(uint8_t));

	copy_sift_image_to_uchar(smooth, image, w + 2, h + 2);

	for (ysrc = s->up, ydest = 0; ysrc < s->down;
			ysrc += sdata->ystep, ydest += sdata->ystep) {
		for (xsrc = s->left, xdest = 0; xsrc < s->right;
				xsrc += sdata->xstep, xdest += sdata->xstep) {
			psrc = xsrc + ysrc * sdata->width;
			if (sdata->umask[psrc] > 0) {
				pdest = xdest + ydest * w;
				psmooth = (xdest + 1) + (ydest + 1) * (w + 2);
				/*
				sobel_win(edge + pdest, w, frame + psrc,
						sdata->width, sdata->xstep,
						sdata->ystep, 192);
						*/
				sobel_win(edge + pdest, w, smooth + psmooth,
						w + 2, sdata->xstep,
						sdata->ystep, 100);
			}
		}
	}

	edge_vertical_profile(&sdata->vprof, edge, w, h, sdata->ystep,
			s->left, s->up);
	window_smooth(&sdata->vprof, 3, 21);
	find_up_down(s, &sdata->vprof, sdata->ystep);

	/*
	save_slope(frame, sdata->width, sdata->height, s->left, s->up,
			s->right, s->down, cnt);
	write_uchar_to_pgm(edge, w, h, "/tmp/%05d-edge.pgm", cnt);
	*/

	cnt++;

	delete smooth;
	delete edge;
	DisallocMatrix(image->pixels, image->rows, image->cols, PERM_POOL);
	free(image);

	return 1;
}

static int find_slope_up_and_down(struct sub_slope *s, struct sub_data *sdata,
		int ystep)
{
	struct profile *prof;
	vector<struct sub_slope> slopes;
	int center;

	prof = &sdata->vprof;

	compute_vertical_profile(sdata, s->left, s->right, s->up, s->down);
	window_smooth(prof, 3, 21);

	slopes = findpeaks(prof, prof->size);
	if (slopes.size() == 0) {
		return 0;
	}
	if (slopes.size() > 0) {
		/*
		center = ystep * slopes.at(0).center;
		s->up = center - 120;
		s->down = center + 120;
		if (s->up < sdata->height / 2) {
			s->up = sdata->height / 2;
		}
		if (s->down > sdata->height - 12) {
			s->down = sdata->height - 12;
		}
		*/
		s->up = ystep * slopes.at(0).left;
		s->down = ystep * slopes.at(0).right;
	}
	return slopes.size();
}

vector<struct sub_slope> compute_region_of_interest(
		struct sub_data *sdata, uint8_t *frame,
		vector<struct sub_slope> prev_slopes)
{
	struct profile *prof;
	vector<struct sub_slope> slopes;
	int slope_width;
	int width;
	int height;
	int xstep;
	int ystep;
	int y, i;
	int ret;

	width = sdata->width;
	height = sdata->height;
	xstep = sdata->xstep;
	ystep = sdata->ystep;
	prof = &sdata->hprof;

	/* Compute the horizontal profile by the sum of vertical umask. */
	compute_horizontal_profile(sdata);

	/* Filter the profile curve for a better detection of peaks. */
	window_smooth(prof, 3, (67 / SUBSAMPLE_STEP_Y) + 1);

	/* Create the slopes based on the peaks found in the profile curve. */
	slopes = findpeaks(prof, prof->size);

	/* Filter slopes with width less than X */
	for (i = 0; i < slopes.size();) {
		slope_width = xstep * (slopes.at(i).right - slopes.at(i).left);
		if (slope_width < 240) {
			/* Remove this slope. */
			slopes.erase(slopes.begin() + i);
			continue;
		}
		i++;
	}

	/* Split the slopes according to the lane's boundary. */
	split_slopes(slopes, sdata->xlane_max, sdata->xlane_min);
	for (i = 0; i < slopes.size();) {
		/* Set the limits of the slope. */

#if 1
		slopes.at(i).up = INT_MAX;
		slopes.at(i).down = INT_MIN;
		for (y = slopes.at(i).left; y < slopes.at(i).right; y++) {
			if (prof->max[y] > slopes.at(i).down) {
				slopes.at(i).down = prof->max[y];
			}
			if (prof->min[y] < slopes.at(i).up) {
				slopes.at(i).up = prof->min[y];
			}
                        /*if (prof->profile[y] < 5) {
                            slopes.at(i).down = -999;
                        }*/
		}
                int max = 0;
		for (y = slopes.at(i).left; y < slopes.at(i).right; y++) {
                    if (prof->profile[y] > max) {
                       max = prof->profile[y];
                    }
                }
                if (max < 30) {
                   slopes.at(i).left = 9999;
                }
#endif

#if 0
		slopes.at(i).up = sdata->wevent.ymin;
		slopes.at(i).down = sdata->wevent.ymax;
#endif

		slopes.at(i).left = xstep * slopes.at(i).left;
		slopes.at(i).center = xstep * slopes.at(i).center;
		slopes.at(i).right = xstep * slopes.at(i).right;
		slopes.at(i).id = 0;
		slopes.at(i).ntries = 0;
		fit_slope_inside_wevent(&slopes.at(i), &sdata->wevent);
		if ((slopes.at(i).up >= slopes.at(i).down) ||
				(slopes.at(i).left >= slopes.at(i).right)) {
			slopes.erase(slopes.begin() + i);
			continue;
		}

		/* Compute the limits up and down for this slope. */
#if 0
		/*
		 * ret = find_slope_up_and_down(&slopes.at(i), sdata, ystep);
		 */
		ret = test_slope_up_and_down(&slopes.at(i), frame, sdata);
		if (0 == ret) {
			slopes.erase(slopes.begin() + i);
			continue;
		}

		/*
		plot_vertical_profile(frame, width, height, sdata->ystep,
				slopes.at(i).left, &sdata->vprof, true);
		*/
#endif

		/* Verify the current ID. */
		check_matching_id(slopes.at(i), prev_slopes);
		if (0 == slopes.at(i).id) {
			sdata->slope_id_cnt++;
			slopes.at(i).id = sdata->slope_id_cnt;
		}
		i++;
	}

	return slopes;
}

void sub_paint_debug(struct sub_data *sdata, uint8_t *frame,
		vector<struct sub_slope> slopes)
{
	const uint8_t vcb[4] = {48, 200, 48, 200};
	const uint8_t vcr[4] = {48, 48, 200, 200};
	struct profile *prof;
	int width;
	int height;
	int xstep;
	int y, i;

	width = sdata->width;
	height = sdata->height;
	xstep = sdata->xstep;
	prof = &sdata->hprof;

	for (i = 0; i < slopes.size(); i++) {
		if (!slopes.at(i).search_plate) {
			continue;
		}
		/*
		mean_chrom_rect(frame, width, height, slopes.at(i).left,
				slopes.at(i).up, slopes.at(i).right,
				slopes.at(i).down, vcb[i % 4], vcr[i % 4]);
				*/
		mean_chrom_rect(frame, width, height, slopes.at(i).left,
				slopes.at(i).up, slopes.at(i).right,
				slopes.at(i).down, 56, 56);
		paint_lum_rect(frame, width, height, slopes.at(i).left,
				slopes.at(i).up, slopes.at(i).right,
				slopes.at(i).down, 5);
		draw_text(frame, width, slopes.at(i).left + 12,
				slopes.at(i).up + 12, 5, 255, 0, " ROI ");
	}

	//sub_plot_profile(frame, width, height, xstep, prof, true);
}

/**
 * Updates the search_plate flag inside each slope.
 * This flag is set when the slope has at least 'slope_min_height' of height
 * and the down value is lower than the 'margin'.  Also, before set it, we
 * check if there is avlicense plate being tracked inside the slope's lane.
 * If not, set the flag slopes.at(k).search_plate.
 */
void update_slope_search_plate(vector<struct sub_slope> &slopes,
		struct lane_limits *lanes, int margin)
{
	const int slope_min_height = 380;
	int height;
	int lane;
	int x;
	int y;
	int k;

	for (k = 0; k < slopes.size(); k++) {
		height = slopes.at(k).down - slopes.at(k).up;
		slopes.at(k).search_plate = false;

		if ((slopes.at(k).down < margin) &&
				(height >= slope_min_height)) {

			x = (slopes.at(k).left + slopes.at(k).right) / 2;
			y = (slopes.at(k).up + slopes.at(k).down) / 2;
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
			lane = lanes->lmask[x + y * lanes->width];
			if ((lane < 1) || (lane > GTRUTH_MAX_LANES)) {
				continue;
			}
			lane--;
			if (!lanes->area[lane].tracking_plate) {
				slopes.at(k).search_plate = true;
			}
		}
	}
}

