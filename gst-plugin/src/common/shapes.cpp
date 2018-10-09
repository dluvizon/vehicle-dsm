#include <stdlib.h>

#include "shapes.h"

#include "config.h"
#include "debug.h"

int compute_lines_polygon_4s(struct polygon_4s *poly)
{
	struct point_double p1;
	struct point_double p2;
	struct point_double p3;
	struct point_double p4;
	struct line_func *line;
	double dx;
	double dy;

	p1 = poly->p[0];
	p2 = poly->p[1];
	p3 = poly->p[2];
	p4 = poly->p[3];

	/* check if the points are in the correct order */
	if (p1.x >= p2.x) {
		print_err("%s: p2.x must be greather than p1.x", __func__);
		return -1;
	}
	if (p1.y >= p3.y) {
		print_err("%s: p3.y must be greather than p1.y", __func__);
		return -1;
	}
	if (p3.x >= p4.x) {
		print_err("%s: p4.x must be greather than p3.x", __func__);
		return -1;
	}
	if (p2.y >= p4.y) {
		print_err("%s: p4.y must be greather than p2.y", __func__);
		return -1;
	}

	/* computing the top line */
	line = &poly->top_refx;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	line->slope = dy / dx;
	line->offset = p1.y - p1.x * line->slope;

	/* computing the bottom line */
	line = &poly->bottom_refx;
	dx = p4.x - p3.x;
	dy = p4.y - p3.y;
	line->slope = dy / dx;
	line->offset = p3.y - p3.x * line->slope;

	/* computing the left line */
	line = &poly->left_refy;
	dy = p3.y - p1.y;
	dx = p3.x - p1.x;
	line->slope = dx / dy;
	line->offset = p1.x - p1.y * line->slope;

	/* computing the right line */
	line = &poly->right_refy;
	dy = p4.y - p2.y;
	dx = p4.x - p2.x;
	line->slope = dx / dy;
	line->offset = p2.x - p2.y * line->slope;

	return 0;
}

void draw_line_to_uchar(uint8_t *img, int width, struct seg_line_int *line)
{
	int x1, x2;
	int y1, y2;
	int dx, dy;
	int i;

	if ((line->init.x == line->end.x) && (line->init.y == line->end.y)) {
		/* line with zero length */
		return;
	}

	if (abs(line->init.x - line->end.x) > abs(line->init.y - line->end.y)) {
		/* this line is more horizontally */
		if (line->init.x > line->end.x) {
			x1 = line->end.x;
			x2 = line->init.x;
			y1 = line->end.y;
			y2 = line->init.y;
		} else {
			x1 = line->init.x;
			x2 = line->end.x;
			y1 = line->init.y;
			y2 = line->end.y;
		}
		dx = x2 - x1;
		dy = y2 - y1;
		for (i = x1; i <= x2; i++) {
			int y = y1 + (i - x1) * dy / dx;
			img[y * width + i] = 255;
		}
	} else {
		/* this line is more vertically */
		if (line->init.y > line->end.y) {
			x1 = line->end.x;
			x2 = line->init.x;
			y1 = line->end.y;
			y2 = line->init.y;
		} else {
			x1 = line->init.x;
			x2 = line->end.x;
			y1 = line->init.y;
			y2 = line->end.y;
		}
		dx = x2 - x1;
		dy = y2 - y1;
		for (i = y1; i <= y2; i++) {
			int x = x1 + (i - y1) * dx / dy;
			img[i * width + x] = 255;
		}
	}
}

void draw_rectanglem_to_uchar(uint8_t *img, int width, struct rectangle_m *rec)
{
	uint8_t *p1;
	uint8_t *p2;
	int i;

	/* first draw the two horizontal lines */
	p1 = &img[rec->ymin * width + rec->xmin];
	p2 = &img[rec->ymax * width + rec->xmin];
	for (i = rec->xmin; i <= rec->xmax; i++, p1++, p2++) {
		*p1 = 127;
		*p2 = 127;
	}
	/* now draw the two vertical lines */
	p1 = &img[rec->ymin * width + rec->xmin];
	p2 = &img[rec->ymin * width + rec->xmax];
	for (i = rec->ymin; i <= rec->ymax; i++, p1 += width, p2 += width) {
		*p1 = 127;
		*p2 = 127;
	}
}

void draw_polygon_to_uchar(uint8_t *img, int width, int height,
		struct polygon_4s *poly)
{
	/* TODO: implemente this function */
}

/**
 * Return 1 if the point is above the line, -1 if is bellow or 0 if the point
 * is in the line.
 */
static int point_to_line_in_x(struct line_func *xline, double x, double y)
{
	double ny;

	ny = x * xline->slope + xline->offset;
	if (y < ny) {
		return 1;
	} else if (y > ny) {
		return -1;
	} else {
		return 0;
	}
}

/**
 * Return 1 if the point is on the left of the line, -1 if is on the right or
 * 0 if the point is in the line.
 */
static int point_to_line_in_y(struct line_func *yline, double x, double y)
{
	double nx;

	nx = y * yline->slope + yline->offset;
	if (x < nx) {
		return 1;
	} else if (x > nx) {
		return -1;
	} else {
		return 0;
	}
}

bool is_point_inside_polygon(struct polygon_4s *poly, double x, double y)
{
	int ret;

	ret = point_to_line_in_x(&poly->bottom_refx, x, y);
	if (-1 == ret) {
		return false;
	}
	ret = point_to_line_in_x(&poly->top_refx, x, y);
	if (1 == ret) {
		return false;
	}
	ret = point_to_line_in_y(&poly->left_refy, x, y);
	if (1 == ret) {
		return false;
	}
	ret = point_to_line_in_y(&poly->right_refy, x, y);
	if (-1 == ret) {
		return false;
	}
	return true;
}

void build_mask_from_polygon(uint8_t *mask, int width, int height,
		struct polygon_4s *poly, uint8_t val)
{
	int p;
	int w;
	int h;

	for (h = 0; h < height; h++) {
		for (w = 0; w < width; w++) {
			if (is_point_inside_polygon(poly, w, h)) {
				p = h * width + w;
				mask[p] = val;
			}
		}
	}
}

