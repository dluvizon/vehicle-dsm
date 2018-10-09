/**
 * @file lanes.h
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 03/02/2015
 */

#ifndef __LANES_H
#define __LANES_H

#include <stdbool.h>

#include "shapes.h"
#include "config.h"

struct lane_area {
	struct polygon_4s poly;
	struct polygon_4s measure;
	bool tracking_plate;
};

struct lane_limits {
	uint8_t *lmask;		/**< Lane area mask. */
	uint8_t *mmask;		/**< Measurement area mask. */
	int width;		/**< Image width. */
	int height;		/**< Image height. */
	struct lane_area area[GTRUTH_MAX_LANES];
};

#endif /* __LANES_H */
