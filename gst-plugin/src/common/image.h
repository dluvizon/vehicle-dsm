#ifndef __IMAGE_H
#define __IMAGE_H

#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>

#include "vector.h"
#include "grouping/grouping.h"

#define INTMAX +9999999
#define INTMIN -9999999

void get_image_dimensions (char *filename, int *nrows, int *ncols);

unsigned char* read_pgm (char *filename, int *nrows, int *ncols);

/**
 * Write a unsigned char image to a pgm file.
 * @param dimg Pointer to an unsigned char image.
 * @param w Image width.
 * @param h Image height.
 * @param fname Filename string folowed by varargs.
 */
void write_uchar_to_pgm(uint8_t *uimg, int w, int h, const char *fname, ...);

void write_nv12_to_jpeg(uint8_t *img, int w, int h, const char *fname, ...);

void write_uchar_to_jpeg(char *destfile, unsigned char *img,
		int width, int height);

void write_uchar_to_ppm (char *filename, unsigned char *img, int nrows, int ncols);

void write_argb_to_ppm (char *filename, unsigned char *img, int nrows, int ncols);

unsigned char* convert_argb_to_gray (unsigned char *image, int nrows, int ncols);

unsigned char* convert_rgb_to_gray (unsigned char *image, int nrows, int ncols);

/**
 * Write a double image to a pgm file.
 * @param dimg Pointer to a double image.
 * @param w Image width.
 * @param h Image height.
 * @param fname Filename string folowed by varargs.
 */
void write_double_to_pgm(double *dimg, int w, int h, const char *fname, ...);

void write_int_to_pgm (char *filename, int *img, int nrows, int ncols);

unsigned char* copy_uchar (unsigned char *src, int nrows, int ncols);

void draw_rectangle (unsigned char *image, int nrows, int ncols, int x, int y, int w, int h, int color);

/**
 * Copy a double image to an unsigned 8bits image and adjust the max and min
 * for to double image to be stored from 0 to 255.
 * @param uimg Pointer to an unsigned char image.
 * @param dimg Pointer to a double image.
 * @param w Image width.
 * @param h Image height.
 */
void image_copy_double_to_uchar(uint8_t *uimg, double *dimg, int w, int h);

void image_overwrite_double_to_uchar(uint8_t *uimg, double *dimg, int w, int h,
		uint8_t threshold);

void copy_uchar_to_double (unsigned char *src, int nrows, int ncols, double *dst);

void fill_double_image (double *frame, int nrows, int ncols, double val);

void invert_image(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *src, unsigned char *dst,
		int nrows, int ncols);

unsigned char *resize_gray_uchar_image_bilinear (unsigned char *pixels, int w, int h, int w2, int h2);

void image_change_interval (double *src, int nrows, int ncols, double vmax, double vmin, int min, int max, unsigned char *dst);

double *resize_gray_uchar_double_image_bilinear (unsigned char *pixels, int w, int h, int h2, int w_mod, int *W2);

unsigned char *resize_color_uchar_image_bilinear (unsigned char *pixels, int w, int h, int w2, int h2, int alpha);

void write_chains_to_uchar (char *filename, unsigned char *img, int nrows, int ncols, Vector** chains, int nsets);

void print_pixel_lum(uint8_t *lum, int w, int h, int x, int y,
		uint8_t gray);

void print_pixel_nv12(uint8_t *chr, int w, int h, int x, int y,
		uint8_t cb, uint8_t cr);

void paint_chrom_rect(uint8_t *frame, int w, int h, int x, int y,
		int maxx, int maxy, int cb, int cr);

void mean_pixel_nv12(uint8_t *chr, int w, int h, int x, int y,
		uint8_t cb, uint8_t cr);

void mean_chrom_rect(uint8_t *frame, int w, int h, int x, int y,
		int maxx, int maxy, uint8_t cb, uint8_t cr);

void paint_lum_rect(uint8_t *frame, int w, int h, int x, int y,
		int maxx, int maxy, int thickness);

#endif /* __IMAGE_H */
