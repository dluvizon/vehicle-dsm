/**
 * file: backgroundsub.h
 */

#ifndef __BACKGROUNDSUB_H
#define __BACKGROUNDSUB_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <list>

#include <stdint.h>

#include "image.h"
#include "shapes.h"
#include "config.h"
#include "lanes.h"

#define SUBSAMPLE_STEP_X 4
#define SUBSAMPLE_STEP_Y 4
#define SUBSAMPLE_AREA (SUBSAMPLE_STEP_X * SUBSAMPLE_STEP_Y)

#if (SUBSAMPLE_STEP_Y < 1) || (SUBSAMPLE_STEP_Y > 67)
# error SUBSAMPLE_STEP_Y must be between 1 to 67
#endif

#define PROF_DELTA_PERC 10	/* Percent */

#define SUB_MHI_DURATION 0.01f	/* Time in milliseconds. */
#define SUB_PERCENTAGE 0.01f
#define SUB_THRESHOLD 44	/* From 0 to 255. */
#define SUB_NFRAMES 2

#ifdef __cplusplus
using namespace std;
#endif

/** @struct sub_slope */
struct sub_slope {
	int mag;	/**< Peak magnitude. */
	int up;		/**< Peak up. */
	int down;	/**< Peak down. */
	int left;	/**< Peak left. */
	int center;	/**< Peak center. */
	int right;	/**< Peak right. */
	unsigned long id;
	bool search_plate;
        int ntries;
};

struct profile {
	int *profile;	/**< Filtered prifile. */
	int *accum;	/**< Accumulator. */
	int *min;
	int *max;
	int *smooth;	/**< Temporary array used to smooth the profile. */
	int *rising;	/**< Guarda a posição de onde há uma borda de subida. */
	int *falling;	/**< Guarda a posição de onde há uma borda de descida.*/
	int size;	/**< Size of profile, min and max. */
};

/** Timeout in frames to match a slope with the groundtruth. */
#define SUB_SLOPE_TIMEOUT_RST 5

/**
 * Statistical data from background subtraction module.
 * This structure shoud be used for one single lane.
 */
struct sub_stat_data {
	FILE *logfile;		/* Log file. */
	unsigned int cnt_flush;

	/* Static variables. */
	unsigned long cnt_slopes;
	unsigned long cnt_vehicles;
	unsigned int slope_timeout;
	float last_speed;
};

/** @struct sub_data */
struct sub_data {
	char *logdir_path;	/**< Path to the log directory of the
				 * background subtraction module. */
	long iframe_last;	/**< Index of the last frame. */
	uint8_t *grid;
	double *mhi;
	uint8_t *umask;
	double **buffers;
	int width;
	int height;
	int threshold;
	struct rectangle_m win;		/**< Window to be computed. */
	struct rectangle_m wevent;	/**< Window of events. */
	int xstep;
	int ystep;
	double area;	/**< Area of the main window, in valid pixels. */
	int iteration;

	struct profile hprof;
	struct profile vprof;
	int xlane_max[3];
	int xlane_min[3];

	/* statistical data */
	struct sub_stat_data stat[GTRUTH_MAX_LANES];
	unsigned long slope_id_cnt;
};

/**
 * Creates a newly allocated struct sub_data.
 * @param width Image width.
 * @param height Image height.
 * @param window Window in the image to process the background.
 * @param xstep Step in X axis.
 * @param ystep Step in Y axis.
 * @return A new struct sub_data.
 */
struct sub_data *background_sub_new(int width, int height,
		struct rectangle_m window, int xstep, int ystep);

/**
 * Create a log directory if it does not exist.
 * @param sdata Pointer to a struct sub_data.
 * @param path String with the path.
 */
void sub_create_logdir(struct sub_data *sdata, const char *path);

/**
 * Write statistical data to log file.
 */
void sub_write_stat_log(struct sub_stat_data *stat);

/**
 * @return A boolean value to indicate changes in background (true for changes,
 * false for no changes).
 */
bool background_subtraction(struct sub_data *sdata, uint8_t *frame, int iframe,
		double timestamp, uint8_t *fdebug);

/**
 * TODO:
 */
vector<struct sub_slope> compute_region_of_interest(
		struct sub_data *sdata, uint8_t *frame,
		vector<struct sub_slope> prev_slopes);

void sub_paint_debug(struct sub_data *sdata, uint8_t *frame,
		vector<struct sub_slope> slopes);

void update_slope_search_plate(vector<struct sub_slope> &slopes,
		struct lane_limits *lanes, int margin);

void debug_save_slope_on_float(uint8_t *frame, int width, int height,
		vector<struct sub_slope> &slopes, int margin);

#endif /* __BACKGROUNDSUB_H */

