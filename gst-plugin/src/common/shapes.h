#ifndef __SHAPES_H
#define __SHAPES_H

#include <stdint.h>

#include "image.h"

struct point_int {
	int x;
	int y;
};

struct point_double {
	double x;
	double y;
};

struct seg_line_int {
	struct point_int init;
	struct point_int end;
};

/** Rectangle based on maximum and minimum. */
struct rectangle_m {
	int xmin;
	int ymin;
	int xmax;
	int ymax;
};

/** Rectangle based in one point reference point and the width and height. */
struct rectangle_p {
	struct point_int p;
	int width;
	int height;
};

struct line_func {
	double offset;
	double slope;
};

/**
 * Four side polygon. This is the simplest case:
 *   +------
 *   |      -----+
 *   |           |
 *    |          |
 *    |          |
 *    +---------+
 * Four vertex, two lines most horizontally and two lines most vertically. 
 */
struct polygon_4s {
	struct point_double p[4];
	struct line_func bottom_refx;
	struct line_func top_refx;
	struct line_func left_refy;
	struct line_func right_refy;
};

int compute_lines_polygon_4s(struct polygon_4s *poly);

void draw_line_to_uchar(uint8_t *img, int width, struct seg_line_int *line);

void draw_rectanglem_to_uchar(uint8_t *img, int width, struct rectangle_m *rec);

void draw_polygon_to_uchar(uint8_t *img, int width, int height,
		struct polygon_4s *poly);

bool is_point_inside_polygon(struct polygon_4s *poly, double x, double y);

void build_mask_from_polygon(uint8_t *mask, int width, int height,
		struct polygon_4s *poly, uint8_t val);

#endif /* __SHAPES_H */
