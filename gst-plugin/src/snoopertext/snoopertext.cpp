/**
 * @file snoopertext.cpp
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 10/11/2014
 */

#define _DLEVEL SNOOPERTEXT_DLEVEL

#include <assert.h>

#include "snoopertext.h"
#include "grouping/grouping.h"
#include "seg/toggle.h"
#include "config.h"
#include "lanes.h"
#include "debug.h"

int snoopertext_load_svm_model(struct snoopertext *s, const char *fname)
{
	assert(s && fname);

	s->model = svm_load_model(fname);
	if (!s->model) {
		return -1;
	}
	if (svm_check_probability_model(s->model) == 0) {
		print_err("Model does not support probabiliy estimates\n");
		return -1;
	}
	return 0;
}

int snoopertext_load_thog_settings(struct snoopertext *s, const char *fname)
{
	assert(s && fname);

	/* load T-HOG settings */
	s->sthog = load_settings(fname);

	return 0;
}

vector<plate> text_filtering(unsigned char *img, int nrows, int ncols,
		Vector** chains, int nsets, struct svm_model *model,
		double *prob, struct_thog sthog, int iframe)
{
	int i, j, k;

	vector<plate> plates;

	for (i = 0; i < nsets; i++) {
		/* Do not show isolated regions: */
		if (vector_count(chains[i]) <= 3) {
			continue;
		}

		int xmin = IMAX;
		int ymin = IMAX;
		int xmax = IMIN;
		int ymax = IMIN;

		/* Computing the chain dimensions based on the regions
		 * inside it: */
		for (j = 0; j < vector_count(chains[i]); j++) {
			region *r = (region *) vector_get(chains[i], j);
			if (r->box[0][0] < xmin) {
				xmin = r->box[0][0];
			}
			if (r->box[1][0] < ymin) {
				ymin = r->box[1][0];
			}
			if (r->box[0][1] > xmax) {
				xmax = r->box[0][1];
			}
			if (r->box[1][1] > ymax) {
				ymax = r->box[1][1];
			}
		}
		int wbox = xmax - xmin;
		int hbox = ymax - ymin;
		if (hbox > wbox) {
			continue;
		}

		if (classify(img, nrows, ncols, xmin, ymin, wbox, hbox,
					model, prob, sthog) > 0) {
			plate p = {xmin, ymin, wbox, hbox};
			p.speed = 0.0;
			p.from = NULL;

			plates.push_back(p);

#if SNOOPERTEXT_DLEVEL > 3
			/* drawing the horizontal line */
			for (k = xmin; k < xmax; k++) {
				img[ncols*ymin + k] = 255;
				img[ncols*ymax + k] = 255;
				img[ncols*(ymin-1) + k] = 0;
				img[ncols*(ymax+1) + k] = 0;
			}
			/* Drawing the vertical line: */
			for (k = ymin; k < ymax; k++) {
				img[ncols*k + xmin] = 255;
				img[ncols*k + xmax] = 255;
				img[ncols*k + xmin - 1] = 0;
				img[ncols*k + xmax + 1] = 0;
			}
			/*
			text_recognition(img, nrows, ncols, xmin, ymin, wbox,
					hbox, iframe, i); 
			*/
#endif
		}
	}

	return plates;
}

/**
 * SnooperText algorithm used to detect license plates.
 * This function returns a vector of license plates.
 */
vector<plate> snoopertext_detection_working(
		unsigned char *image,
		struct snoopertext *s, int nrows, int ncols, int iframe,
		int xstep, vector<struct sub_slope> &slopes,
		int margin_X, int margin_Y)
{
	vector<plate> plates;
	Vector ncc;
	Vector icc;
        Vector** ochains;

	int nsets;
	int isets;
	int i;

     /*Parameter settings: */
	/* Toggle segmentation settings: */
        int mask_size = 19;                    /*change*/
	int contrast = 18;                     /*change*/
	int percentage = 86;                   /*change*/

	/* Geometric filtering settings: */
	double rmin = 3;                       /*change*/
	double rmax = rmin * 5;                /*fixed*/
	double fmin = 0.1;                     /*fixed*/

	/* Character grouping settings: */
        double t1 = 0.7;                       /*fixed*/
        double t2 = 1.1;                       /*fixed*/
        double t3 = 0.4;                       /*fixed*/
     /*End settings*/

	int cter = 0;

	/* Clear all buffers. */
	memset(s->nmin, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->nmax, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->imin, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->imax, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->nseg, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->iseg, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->ntmp, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->itmp, 0, nrows * ncols * sizeof(unsigned char));
	memset(s->invert, 0, nrows * ncols * sizeof(unsigned char));

	for (i = 0; i < slopes.size(); i++) {
		if (!slopes.at(i).search_plate) {
			continue;
		}
		cter++;
	}
	if (cter == 0) {
		goto end;
	}

	for (i = 0; i < slopes.size(); i++) {
		if (!slopes.at(i).search_plate) {
			continue;
		}
                slopes.at(i).left = slopes.at(i).left - margin_X;
                slopes.at(i).right = slopes.at(i).right + margin_X;
                slopes.at(i).up = slopes.at(i).up - margin_Y;
                slopes.at(i).down = slopes.at(i).down + margin_Y;

		/* Compute the inverted image in each slope area. */
		invert_image(slopes.at(i).left, slopes.at(i).up,
				slopes.at(i).right, slopes.at(i).down,
				image, s->invert, nrows, ncols);
	}

	/**
	 * Compute images after dilation and erision for both normal and
	 * inverted images.  This will generate the images following images:
	 *   s->nmax: Normal maximum (dilate)
	 *   s->imax: Inverted maximum (dilate
	 *   s->nmin: Normal minimum (erode)
	 *   s->imin: Inverted minimum (erode)
	 */
	for (i = 0; i < slopes.size(); i++) {
		if (!slopes.at(i).search_plate) {
			continue;
		}
		dilate(slopes.at(i).left, slopes.at(i).up,
				slopes.at(i).right, slopes.at(i).down,
				image, s->ntmp, s->invert, s->itmp,
				ncols, nrows, mask_size, HORZ);
		dilate(slopes.at(i).left, slopes.at(i).up,
				slopes.at(i).right, slopes.at(i).down,
				s->ntmp, s->nmax, s->itmp, s->imax,
				ncols, nrows, mask_size, VERT);
		erode(slopes.at(i).left, slopes.at(i).up,
				slopes.at(i).right, slopes.at(i).down,
				image, s->ntmp, s->invert, s->itmp,
				ncols, nrows, mask_size, HORZ);
		erode(slopes.at(i).left, slopes.at(i).up,
				slopes.at(i).right, slopes.at(i).down,
				s->ntmp, s->nmin, s->itmp, s->imin,
				ncols, nrows, mask_size, VERT);
	}

	/**
	 * Toggle function: implement the segmentation step.
	 * Output: buffers s->nseg and s->iseg.
	 */
	for (i = 0; i < slopes.size(); i++) {
		if (!slopes.at(i).search_plate) {
			continue;
		}
		toggle(slopes.at(i).left, slopes.at(i).up,
				slopes.at(i).right, slopes.at(i).down,
				image, s->nmin, s->nmax, s->nseg,
				s->invert, s->imin, s->imax, s->iseg,
				ncols, nrows, mask_size, contrast, percentage);
	}

#if 0
	/* If enabled, write the images resulted from segmentation. */
	write_uchar_to_pgm(s->nseg, ncols, nrows, "/tmp/snooper/nseg_%05d.pgm",
			iframe);
	write_uchar_to_pgm(s->iseg, ncols, nrows, "/tmp/snooper/iseg_%05d.pgm",
			iframe);
#endif

	ncc = find_connected_components_working(s->nseg, nrows, ncols,
			rmin, rmax, fmin, 0, 0, ncols, nrows);
	icc = find_connected_components_working(s->iseg, nrows, ncols,
			rmin, rmax, fmin, 0, 0, ncols, nrows);

        ochains = make_chains_and_merge(ncc, icc, t1, t2, t3, &nsets);
	plates = text_filtering(image, nrows, ncols, ochains, nsets + isets,
			s->model, s->prob, s->sthog, iframe);

	chain_free(ochains, nsets + isets);

end:
	return plates;
}

