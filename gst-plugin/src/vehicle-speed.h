/**
 * @file vehicle-speed.h
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 29/01/2014
 */

#ifndef __VEHICLE_SPEED_H
#define __VEHICLE_SPEED_H

#include <stdio.h>
#include <stdbool.h>

#include "klt/klt.h"

#include "backgroundsub/backgroundsub.h"
#include "gtruth/vehicles.h"
#include "snoopertext/snoopertext.h"
#include "lanes.h"

/* Input files. */
#define SVM_MODEL_FILE "./input/model.svm"
#define THOG_SETTINGS_FILE "./input/1_7_9.txt"

struct dsm_vs_data {
	int last;

	unsigned long object_counter;
	const char *out_path;

	/* Image parameters. */
	int img_width;
	int img_height;

	/* KLT */
	KLT_TrackingContext tc;
	KLT_TrackingContext tc_gt;
	vector<KLT_FeatureList> features;
	vector<KLT_FeatureList> features_gt;
	struct klt_buffers kltbufs;

	struct lane_limits lanes;

	struct sub_data *background_sub;
	struct vehicle_table *vt;
	struct snoopertext snooper;

	bool should_save_frame;
};

#ifndef max
# define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
# define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#endif /* __VEHICLE_SPEED_H */

