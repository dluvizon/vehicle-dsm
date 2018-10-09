/**
 * @file snoopertext.h
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 10/11/2014
 */

#ifndef __SNOOPERTEXT_H
#define __SNOOPERTEXT_H

#include <vector>

#include "backgroundsub/backgroundsub.h"
#include "seg/toggle.h"
#include "svm/svm.h"
#include "thog/thog.h"
#include "image.h"

struct snoopertext {
	struct svm_model *model;
	struct_thog sthog;
	unsigned char *nmin;
	unsigned char *nmax;
	unsigned char *imin;
	unsigned char *imax;
	unsigned char *nseg;
	unsigned char *iseg;
	unsigned char *ntmp;
	unsigned char *itmp;
	unsigned char *invert;
	unsigned char *edges;
	unsigned char *dilate;
	unsigned char *enhanced;
	double *prob;
};

int snoopertext_load_svm_model(struct snoopertext *s, const char *fname);

int snoopertext_load_thog_settings(struct snoopertext *s, const char *fname);

vector<plate> text_filtering(unsigned char *img, int nrows, int ncols,
		Vector** chains, int nsets, struct svm_model *model,
		double *prob, struct_thog sthog, int iframe);

vector<plate> snoopertext_detection_working(
		unsigned char *image,
		struct snoopertext *s,
		int nrows,
		int ncols,
		int iframe,
		int xstep,
		vector<struct sub_slope> &slopes,
		int margin_X,
		int margin_Y);

#endif /* __SNOOPERTEXT_H */
