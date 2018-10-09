#include <assert.h>
#include <stdio.h>
#include <float.h>

#include <jpeglib.h>

#include "image.h"
#include "debug.h"

/* */
void get_image_dimensions (char *filename, int *nrows, int *ncols)
{
   FILE *file = fopen (filename, "r");

   char *type = (char *)malloc(3*sizeof(char));

   if (!fscanf (file, "%s", type)) {
	   return;
   }

   free(type);

   if (!fscanf (file, "%d %d", ncols, nrows)) {
	   return;
   }

   fclose (file);
}   

/**
 * FIXME: This function have bugs, it is not working!
 * Maybe something in scanf() call. It is not checking the maximum value
 * of the grayscale pixels when reading file.
 */
unsigned char* read_pgm (char *filename, int *nrows, int *ncols) {

   FILE *file = fopen (filename, "r");

   char *type = (char *)malloc(3*sizeof(char));

   if (!fscanf (file, "%s", type)) {
	   return NULL;
   }

   free(type);

   if (!fscanf (file, "%d %d", ncols, nrows)) {
	   return NULL;
   }

   int size = (*nrows) * (*ncols);

   unsigned char *img = (unsigned char *)malloc(size * sizeof(char));

   int i;
   for (i = 0; i < size; i++) {
      int pixel;
      if (fscanf (file, "%d", &pixel)) {
	      return NULL;
      }
      img[i] = (unsigned char)pixel;
   }
   fclose (file);
   return img;
}

void write_uchar_to_pgm(uint8_t *uimg, int w, int h, const char *fname, ...)
{
	char filename[256];
	va_list argptr;
	FILE *file;
	int size;
	int i;

	va_start(argptr, fname);
	vsnprintf(filename, sizeof(filename), fname, argptr);
	va_end(argptr);

	file = fopen(filename, "w");
	if (!file) {
		print_err("could not open file %s", filename);
		return;
	}

	fprintf(file, "P2\n");
	fprintf(file, "%d %d\n", w, h);
	fprintf(file, "255\n");
	size = w * h;
	for (i = 0; i < size; i++) {
		fprintf(file, "%d ", (int) uimg[i]);
		if (i % 12 == 0) {
			fprintf(file, "\n");
		}
	}

	fclose(file);
}

static void image_convert_nv12_to_420p(uint8_t *cb, uint8_t *cr,
		uint8_t *nv12, unsigned int ysize)
{
	unsigned int ch_size;
	uint8_t *nv12_end;

	ch_size = ysize / 2;

	nv12 = nv12 + ysize;
	nv12_end = nv12 + ch_size;

	while (nv12 < nv12_end - 1) {
		*cb++ = *nv12++;
		*cr++ = *nv12++;
	}
}

static void jpeg_encode_color_nv12(const char *destfile, uint8_t *imgp,
		int w, int h)
{
	static uint8_t buffcb[1920 * 1080 / 4];
	static uint8_t buffcr[1920 * 1080 / 4];
	FILE *file;
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPARRAY data[3];
	JSAMPROW y[2 * DCTSIZE];
	JSAMPROW cb[DCTSIZE];
	JSAMPROW cr[DCTSIZE];
	uint8_t *img_y;
	uint8_t *img_cb;
	uint8_t *img_cr;
	unsigned int offset;
	unsigned int i;

	file = fopen(destfile, "w");
	if (!file) {
		fprintf(stderr, "Could nop open file '%s' for writting\n",
				destfile);
		exit(1);
	}

	image_convert_nv12_to_420p(buffcb, buffcr, imgp, w * h);
	data[0] = y;
	data[1] = cb;
	data[2] = cr;

	/* create a jpeg compress and configure it */
	jpeg_create_compress(&cinfo);
	cinfo.err = jpeg_std_error(&jerr);
	cinfo.image_width = w;
	cinfo.image_height = h;
	/* configure to look for luminance, cb and cr */
	cinfo.input_components = 3;
	jpeg_set_defaults(&cinfo);
	jpeg_set_colorspace(&cinfo, JCS_YCbCr);

	/* sample factor according to fourcc nv12 */
	cinfo.comp_info[0].h_samp_factor = 2;	/* y */
	cinfo.comp_info[0].v_samp_factor = 2;
	cinfo.comp_info[1].h_samp_factor = 1;	/* cb */
	cinfo.comp_info[1].v_samp_factor = 1;
	cinfo.comp_info[2].h_samp_factor = 1;	/* cr */
	cinfo.comp_info[2].v_samp_factor = 1;

#if JPEG_LIB_VERSION >= 70
	cinfo.do_fancy_downsampling = FALSE;
	cinfo.dct_method = JDCT_FASTEST;
	cinfo.smoothing_factor = 0;
#endif
	jpeg_set_quality(&cinfo, 90, TRUE);
	cinfo.raw_data_in = TRUE;

	/* configure the jpeg destination */
	jpeg_stdio_dest(&cinfo, file);
	jpeg_start_compress(&cinfo, TRUE);

	/* encoding */
	while (cinfo.next_scanline < cinfo.image_height) {
		img_y = imgp + (cinfo.next_scanline * w);
		offset = (cinfo.next_scanline * w) / 4;
		img_cb = buffcb + offset;
		img_cr = buffcr + offset;

		for (i = 0; i < 2 * DCTSIZE; i++) {
			y[i] = img_y + (i * w);
			if (i % 2 == 0) {
				offset = (i * w / 4);
				cb[i / 2] = img_cb + offset;
				cr[i / 2] = img_cr + offset;
			}
		}
		jpeg_write_raw_data(&cinfo, data, 2 * DCTSIZE);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	fclose(file);
}

void write_nv12_to_jpeg(uint8_t *img, int w, int h, const char *fname, ...)
{
	char filename[256];
	va_list argptr;
	FILE *file;
	int size;
	int i;

	va_start(argptr, fname);
	vsnprintf(filename, sizeof(filename), fname, argptr);
	va_end(argptr);

	jpeg_encode_color_nv12(filename, img, w, h);
}

void write_uchar_to_jpeg(char *destfile, unsigned char *img,
		int width, int height)
{
	FILE *file;
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW y[DCTSIZE];
	int row_stride;

	file = fopen(destfile, "w");
	if (!file) {
		fprintf(stderr, "Could nop open file '%s' for writting\n",
				destfile);
		exit(1);
	}

	jpeg_create_compress(&cinfo);
	cinfo.image_width = width;
	cinfo.image_height = height;
	cinfo.input_components = 1;
	cinfo.in_color_space = JCS_GRAYSCALE;
	row_stride = width;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_set_defaults(&cinfo);
#if JPEG_LIB_VERSION >= 70
	cinfo.do_fancy_downsampling = FALSE;
	cinfo.dct_method = JDCT_FASTEST;
	cinfo.smoothing_factor = 0;
#endif
	jpeg_set_quality(&cinfo, 70, TRUE);

	/* configure the jpeg destination */
	jpeg_stdio_dest(&cinfo, file);
	jpeg_start_compress(&cinfo, TRUE);

	/* encoding */
	while (cinfo.next_scanline < cinfo.image_height) {
		y[0] = img + (cinfo.next_scanline * row_stride);
		jpeg_write_scanlines(&cinfo, y, 1);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	fclose(file);
}

void write_uchar_to_ppm (char *filename, unsigned char *img, int nrows, int ncols) {

   FILE *file = fopen (filename, "w");

   fprintf (file, "P3\n");

   fprintf (file, "%d %d\n", ncols, nrows);

   fprintf (file, "255\n");

   int i, j;
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
         int r = img[3 * i * ncols + 3 * j + 0]; /*red*/
         int g = img[3 * i * ncols + 3 * j + 1]; /*green*/
         int b = img[3 * i * ncols + 3 * j + 2]; /*blue*/
         //int a = img[4 * i * ncols + 4 * j + 3]; /*blue*/
         fprintf (file, "%d %d %d ", r, g, b); 
         if (j % 12 == 0) { fprintf (file, "\n"); }
      }
   }
   fclose (file);
}

/* */
void write_argb_to_ppm (char *filename, unsigned char *img, int nrows, int ncols) {

   FILE *file = fopen (filename, "w");

   fprintf (file, "P3\n");

   fprintf (file, "%d %d\n", ncols, nrows);

   fprintf (file, "255\n");

   int i, j;
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
         int r = img[4 * i * ncols + 4 * j + 0]; /*red*/
         int g = img[4 * i * ncols + 4 * j + 1]; /*green*/
         int b = img[4 * i * ncols + 4 * j + 2]; /*blue*/
	 /*
         int a = img[4 * i * ncols + 4 * j + 3];
	 */
         fprintf (file, "%d %d %d ", r, g, b); 
         if (j % 12 == 0) { fprintf (file, "\n"); }
      }
   }
   fclose (file);
}


void image_change_interval (double *src, int nrows, int ncols, double vmax, double vmin, int min, int max, unsigned char *dst) {

   int i;
   int size = nrows * ncols;

   if (vmax == vmin) {
      vmax = -HUGE_VAL; vmin = +HUGE_VAL;
      for (i = 0; i < size; i++) {
         vmax = fmax (vmax, src[i]);
         vmin = fmin (vmin, src[i]);
      }
   }
     
   for (i = 0; i < size; i++) {
      dst[i] = ( (max - min)*(src[i] - vmin) )/(vmax - vmin) + min;
   } 
}

/* */
void write_int_to_pgm (char *filename, int *img, int nrows, int ncols) {

   FILE *file = fopen (filename, "w");

   fprintf (file, "P2\n");

   fprintf (file, "%d %d\n", ncols, nrows);

   fprintf (file, "255\n");

   int size = nrows * ncols;

   int i;
   for (i = 0; i < size; i++) {
      fprintf (file, "%d ", img[i]); 
      if (i % 12 == 0) { fprintf (file, "\n"); }
   }
   fclose (file);
}

void write_double_to_pgm(double *dimg, int w, int h, const char *fname, ...)
{
	char filename[256];
	va_list argptr;
	FILE *file;
	int size;
	int i;

	va_start(argptr, fname);
	vsnprintf(filename, sizeof(filename), fname, argptr);
	va_end(argptr);

	file = fopen(filename, "w");
	if (!file) {
		print_err("could not open file %s", filename);
		return;
	}

	fprintf(file, "P2\n");
	fprintf(file, "%d %d\n", w, h);
	fprintf(file, "255\n");
	size = w * h;
	for (i = 0; i < size; i++) {
		fprintf(file, "%d ", (int) dimg[i]);
		if (i % 12 == 0) {
			fprintf(file, "\n");
		}
	}

	fclose(file);
}

/* */
unsigned char* convert_argb_to_gray (unsigned char *image, int nrows, int ncols) {

   int x, y;

   int size = nrows * ncols;

   unsigned char *out = (unsigned char *)malloc(size * sizeof(unsigned char));

   for (y = 0; y < nrows; y++) {
      for (x = 0; x < ncols; x++) {
         int r = image[4 * y * ncols + 4 * x + 0]; /*red*/
         int g = image[4 * y * ncols + 4 * x + 1]; /*green*/
         int b = image[4 * y * ncols + 4 * x + 2]; /*blue*/
#if 0
         int a = image[4 * y * ncols + 4 * x + 3]; /*alpha*/
#endif
         //r /= 255.0; g /= 255.0; b /= 255.0;
         unsigned char gray = (unsigned char)(0.299*r + 0.587*g + 0.114*b);
         out[y * ncols + x] = gray;
#if 0
         a = 0; /*The alpha channel is not used.*/
#endif
      }
   }
   return out;
}

/* */
void draw_rectangle (unsigned char *image, int nrows, int ncols, int x, int y, int w, int h, int color) {
   int k;
   for (k = x; k < (x + w); k++) {
      image[y * ncols + k] = color;
      image[(y + h)* ncols + k] = color;
   }
   for (k = y; k < (y + h); k++) {
      image[k * ncols + x] = color;
      image[k * ncols + (x + w)] = color;
   }
}

/* */
unsigned char* convert_rgb_to_gray (unsigned char *image, int nrows, int ncols) {

   int x, y;

   int size = nrows * ncols;

   unsigned char *out = (unsigned char *)malloc(size * sizeof(unsigned char));

   for (y = 0; y < nrows; y++) {
      for (x = 0; x < ncols; x++) {
         int r = image[3 * y * ncols + 3 * x + 0]; /*red*/
         int g = image[3 * y * ncols + 3 * x + 1]; /*green*/
         int b = image[3 * y * ncols + 3 * x + 2]; /*blue*/
         //r /= 255.0; g /= 255.0; b /= 255.0;
         unsigned char gray = (unsigned char)(0.299*r + 0.587*g + 0.114*b);
         out[y * ncols + x] = gray;
      }
   }
   return out;
}

unsigned char *copy_uchar(unsigned char *src, int nrows, int ncols)
{
	int size = nrows * ncols;
	unsigned char *dst;
	int i;

	dst = (unsigned char *) malloc(size * sizeof(unsigned char));
	for (i = 0; i < size; i++) {
		dst[i] = src[i];
	}

	return dst;
}

void image_copy_double_to_uchar(uint8_t *uimg, double *dimg, int w, int h)
{
	double *pd;
	double max;
	double min;
	double coef;
	int size;

	size = w * h;
	max = DBL_MIN;
	min = DBL_MAX;
	pd = dimg;
	while (size--) {
		if (*pd > max) {
			max = *pd;
		}
		if (*pd < min) {
			min = *pd;
		}
		pd++;
	}

	size = w * h;
	coef = 255. / (max - min);
	while (size--) {
		*uimg++ = (uint8_t) ((*dimg++ - min) * coef);
	}
}

void image_overwrite_double_to_uchar(uint8_t *uimg, double *dimg, int w, int h,
		uint8_t threshold)
{
	double *pd;
	double max;
	double min;
	double coef;
	uint8_t pix;
	int size;

	size = w * h;
	max = DBL_MIN;
	min = DBL_MAX;
	pd = dimg;
	while (size--) {
		if (*pd > max) {
			max = *pd;
		}
		if (*pd < min) {
			min = *pd;
		}
		pd++;
	}

	printf("MIN/MAX = %.1f/%.1f\n", min, max);

	size = w * h;
	coef = 255. / (max - min);
	while (size--) {
		pix = (uint8_t) ((*dimg++ - min) * coef);
		if (pix >= threshold) {
			*uimg = 255;
		}
		uimg++;
	}
}

/*Assumes that both images {src} and {dst} have the same size {nrows*cols}. */
void copy_uchar_to_double (unsigned char *src, int nrows, int ncols, double *dst) {

   int size = nrows * ncols;

   int i;
   for (i = 0; i < size; i++) {
      dst[i] = src[i];
   }
}

/* */
void fill_double_image (double *frame, int nrows, int ncols, double val) {

   int size = nrows * ncols;

   int i;
   for (i = 0; i < size; i++) {
      frame[i] = val;
   }
}

/* */
void invert_image(int xmin, int ymin, int xmax, int ymax,
		const unsigned char *src, unsigned char *dst,
		int nrows, int ncols)
{
	int x, y;
	int i;

	for (y = ymin; y < ymax; y++) {
		for (x = xmin; x < xmax; x++) {
			i = y * ncols + x;
			dst[i] = 255 - src[i];
		}
	}
}

unsigned char *resize_gray_uchar_image_bilinear (unsigned char *pixels, int w, int h, int w2, int h2) {

    unsigned char *temp = (unsigned char *)malloc(w2*h2*sizeof(unsigned char));
    int A, B, C, D, x, y, index, gray ;
    float x_ratio = ((float)(w-1))/w2 ;
    float y_ratio = ((float)(h-1))/h2 ;
    float x_diff, y_diff;
    int offset = 0 ;
    int i, j;
    for (i=0;i<h2;i++) {
        for (j=0;j<w2;j++) {
            x = (int)(x_ratio * j) ;
            y = (int)(y_ratio * i) ;
            x_diff = (x_ratio * j) - x ;
            y_diff = (y_ratio * i) - y ;
            index = y*w+x ;

            // range is 0 to 255 thus bitwise AND with 0xff
            A = pixels[index] & 0xff ;
            B = pixels[index+1] & 0xff ;
            C = pixels[index+w] & 0xff ;
            D = pixels[index+w+1] & 0xff ;
            
            // Y = A(1-w)(1-h) + B(w)(1-h) + C(h)(1-w) + Dwh
            gray = (int)(
                    A*(1-x_diff)*(1-y_diff) +  B*(x_diff)*(1-y_diff) +
                    C*(y_diff)*(1-x_diff)   +  D*(x_diff*y_diff)
                    ) ;

            temp[offset++] = gray ;                                   
        }
    }
    return temp ;
}

double *resize_gray_uchar_double_image_bilinear (unsigned char *pixels, int w, int h, int h2, int w_mod, int *W2) {

    double y_ratio = (double)(h2)/(double)(h);
    int w2 = (int)(floor(w*y_ratio/w_mod + 0.5))*w_mod;
    *W2 = w2; 
    assert((w2 % w_mod) == 0);
    double x_ratio = (double)(w2)/(double)(w); 

    double *temp = (double *)malloc(w2*h2*sizeof(double));
    int A, B, C, D, x, y, index, gray ;
    double x_diff, y_diff;
    int offset = 0 ;
    int i, j;
    for (i=0;i<h2;i++) {
        for (j=0;j<w2;j++) {
            x = (int)(x_ratio * j) ;
            y = (int)(y_ratio * i) ;
            x_diff = (x_ratio * j) - x ;
            y_diff = (y_ratio * i) - y ;
            index = y*w+x ;

            // range is 0 to 255 thus bitwise AND with 0xff
            A = pixels[index] & 0xff ;
            B = pixels[index+1] & 0xff ;
            C = pixels[index+w] & 0xff ;
            D = pixels[index+w+1] & 0xff ;
            
            gray =  ( A*(1-x_diff)*(1-y_diff) +  B*(x_diff)*(1-y_diff) + C*(y_diff)*(1-x_diff)   +  D*(x_diff*y_diff) );

            temp[offset++] = gray ;                                   
        }
    }
    return temp ;
}

/**/
unsigned char *resize_color_uchar_image_bilinear (unsigned char *pixels, int w, int h, int w2, int h2, int alpha) {

    int nchannels;

    if (alpha) { nchannels = 4; }

    else { nchannels = 3; }

    int nch = 3;
    unsigned char *temp = (unsigned char *)malloc(nch*w2*h2*sizeof(unsigned char));
    int A, B, C, D, x, y, index, gray ;
    float x_ratio = ((float)(w-1))/w2 ;
    float y_ratio = ((float)(h-1))/h2 ;
    float x_diff, y_diff;
    int offset = 0 ;
    int i, j;
    //printf("lala\n");
    for (i=0;i<h2;i++) {
        for (j=0;j<w2;j++) {
            x = (int)(x_ratio * j);
            y = (int)(y_ratio * i);
            x_diff = (x_ratio * j) - x;
            y_diff = (y_ratio * i) - y;
            int c;
            for (c = 0; c < nch; c++) {
               index = nchannels * y * w + nchannels * x + c;
               //if (index < nchannels*w*h) {
                  //printf("aqui 00\n");
                  A = pixels[nchannels * y * w + nchannels * x + c];
               //}   
               //if (index+1 < nchannels*w*h) {
                  //printf("aqui 01\n");
                  B = pixels[nchannels * y * w + nchannels * (x+1) + c];
               //}   
               //if (index+w < nchannels*w*h) {
                  //printf("aqui 02\n");
                  C = pixels[nchannels * (y+1) * w + nchannels * x + c];
               //}   
               //if (index+w+1 < nchannels*w*h) {
                 // printf("aqui 03\n");
                  D = pixels[nchannels * (y+1) * w + nchannels * (x+1) + c];
               //}  
               //printf("aqui 5\n");
               temp[offset++] = ( A*(1-x_diff)*(1-y_diff) +  B*(x_diff)*(1-y_diff) + C*(y_diff)*(1-x_diff) + D*(x_diff*y_diff));
               //printf("aqui 6\n");
            }
        }
    }
    //printf("aqui 12\n");
    return temp ;
}

/* */
#if 0
void write_chains_to_uchar (char *filename, unsigned char *img, int nrows, int ncols, Vector** chains, int nsets) {

   int i, j, k;

   unsigned char *out = copy_uchar (img, nrows, ncols);

   for (i = 0; i < nsets; i++) {

      /*Do not show isolated regions: */
      if (vector_count(chains[i]) <= 1) { continue; }

      int xmin = INTMAX, ymin = INTMAX, xmax = INTMIN, ymax = INTMIN;

      /*Computing the chain dimensions based on the regions inside it: */
      for (j = 0; j < vector_count(chains[i]); j++) {
         region *r = (region *)vector_get(chains[i], j);
         if (r->box[0][0] < xmin) { xmin = r->box[0][0]; }
         if (r->box[1][0] < ymin) { ymin = r->box[1][0]; }
         if (r->box[0][1] > xmax) { xmax = r->box[0][1]; }
         if (r->box[1][1] > ymax) { ymax = r->box[1][1]; }

      }
      int wbox = xmax - xmin;
      int hbox = ymax - ymin;
      if (hbox > wbox) { continue; }
      /*Drawing the horizontal line: */
      for (k = xmin; k < xmax; k++) {
         out[ncols*ymin + k] = 255;
         out[ncols*ymax + k] = 255;
      }
      /*Drawing the vertical line: */
      for (k = ymin; k < ymax; k++) {
         out[ncols*k + xmin] = 255;
         out[ncols*k + xmax] = 255;
      }
   }
   write_uchar_to_pgm (filename, out, nrows, ncols);
   free(out);
}
#endif

void print_pixel_lum(uint8_t *lum, int w, int h, int x, int y,
		uint8_t gray)
{
	int offset;

	assert((w > x) && (h > y));

	offset = y * w + x;
	lum[offset] = gray;
}

void print_pixel_nv12(uint8_t *chr, int w, int h, int x, int y,
		uint8_t cb, uint8_t cr)
{
	int chy;
	int chx;
	int p;

	assert((w > x) && (h > y));

	chx = (x / 2) * 2;
	chy = y / 2;
	p = chy * w + chx;
	chr[p] = cb;
	chr[p + 1] = cr;
}

void mean_pixel_nv12(uint8_t *chr, int w, int h, int x, int y,
		uint8_t cb, uint8_t cr)
{
	int chy;
	int chx;
	int p;

	assert((w > x) && (h > y));

	chx = (x / 2) * 2;
	chy = y / 2;
	p = chy * w + chx;

	chr[p] = (uint8_t) (((int) chr[p] + (int) cb) / 2);
	chr[p + 1] = (uint8_t) (((int) chr[p + 1] + (int) cr) / 2);
}

void paint_chrom_rect(uint8_t *frame, int w, int h, int x, int y,
		int maxx, int maxy, int cb, int cr)
{
	uint8_t *chrom;
	int i, j;

	if (maxx > w - 5) {
		maxx = w - 5;
	}
	if (maxy > h - 5) {
		maxy = h - 5;
	}
	if (x < 4) {
		x = 4;
	}
	if (y < 4) {
		y = 4;
	}
	
	chrom = frame + w * h;
	for (j = y; j < maxy; j += 2) {
		for (i = x; i < maxx; i += 2) {
			print_pixel_nv12(chrom, w, h, i, j, cb, cr);
		}
	}
}

void mean_chrom_rect(uint8_t *frame, int w, int h, int x, int y,
		int maxx, int maxy, uint8_t cb, uint8_t cr)
{
	uint8_t *chrom;
	int i, j;

	if (maxx > w - 5) {
		maxx = w - 5;
	}
	if (maxy > h - 5) {
		maxy = h - 5;
	}
	if (x < 4) {
		x = 4;
	}
	if (y < 4) {
		y = 4;
	}

	chrom = frame + w * h;
	for (j = y; j < maxy; j += 2) {
		for (i = x; i < maxx; i += 2) {
			mean_pixel_nv12(chrom, w, h, i, j, cb, cr);
		}
	}
}

void paint_lum_rect(uint8_t *frame, int w, int h, int x, int y,
		int maxx, int maxy, int thickness)
{
	int ix;
	int iy;

	/* Horizontal top line. */
	for (iy = y; iy < y + thickness; iy++) {
		for (ix = x; ix < maxx; ix++) {
			print_pixel_lum(frame, w, h, ix, iy, 255);
		}
	}
	/* Horizontal bottom line. */
	for (iy = maxy; iy > maxy - thickness; iy--) {
		for (ix = x; ix < maxx; ix++) {
			print_pixel_lum(frame, w, h, ix, iy, 255);
		}
	}
	/* Vertical left line. */
	for (ix = x; ix < x + thickness; ix++) {
		for (iy = y; iy < maxy; iy++) {
			print_pixel_lum(frame, w, h, ix, iy, 255);
		}
	}
	/* Vertical right line. */
	for (ix = maxx; ix > maxx - thickness; ix--) {
		for (iy = y; iy < maxy; iy++) {
			print_pixel_lum(frame, w, h, ix, iy, 255);
		}
	}
}

