#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <cctype>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <iomanip>

#include "hmatrix.h"

using namespace std;
using namespace cv;

#define PIXELS_PER_M 60.0

struct point {
	float x;
	float y;
};

/* Inverse perspective mapping. */
struct quadrangle {
	struct point p_top_l;
	struct point p_top_r;
	struct point p_bot_l;
	struct point p_bot_r;
	/* Height of the reference rect in the real world, in meters. */
	float ref_height;
	/* Width of the reference rect in the real world, in meters. */
	float ref_width;
	/* Desired number of pixels per meter in the inverse perspective
	 * mapping. */
	float pixels_per_m;
	/* Top position of the reference rect in the inverse perspective
	 * mapped image. */
	float ipm_top;
	/* Left position of the reference rect in the inverse perspective
	 * mapped image. */
	float ipm_left;
};

static void print_hmatrix(struct hmatrix *mat, int num)
{
	cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
	cout.precision(8);
	cout << "{" << endl << "\t" <<
		"/* faixa " << num << ": */" << endl << "\t" <<
		mat->a0 << ", " << mat->a1 << ", " << mat->a2 << "," <<
		endl << "\t" <<
		mat->b0 << ", " << mat->b1 << ", " << mat->b2 << "," <<
		endl << "\t" <<
		mat->c0 << ", " << mat->c1 << ", " << mat->c2 <<
		endl << "},";
}

static void save_hmatrix(FILE *file, struct hmatrix *mat)
{
	fprintf(file, "%.8f\t%.8f\t%.8f\n", mat->a0, mat->a1, mat->a2);
	fprintf(file, "%.8f\t%.8f\t%.8f\n", mat->b0, mat->b1, mat->b2);
	fprintf(file, "%.8f\t%.8f\t%.8f\n", mat->c0, mat->c1, mat->c2);
	fprintf(file, "\n");
}

static Mat generate_hmatrix(struct hmatrix *hmat,
		const struct quadrangle *quad)
{
	Point2f ref_points[4];
	Point2f ipm_points [4];
	Mat ipm_matrix;

	/* Reference points in the input image.
	 * We use these to perform an inverse perspective mapping. */

	ref_points[0] = Point2f(quad->p_top_l.x, quad->p_top_l.y);
	ref_points[1] = Point2f(quad->p_top_r.x, quad->p_top_r.y);
	ref_points[2] = Point2f(quad->p_bot_l.x, quad->p_bot_l.y);
	ref_points[3] = Point2f(quad->p_bot_r.x, quad->p_bot_r.y);

	/* Reference points in the inverse perspective mapping. */
	ipm_points[0] = Point2f(quad->ipm_left, quad->ipm_top);
	ipm_points[1] = Point2f(quad->ipm_left +
			quad->pixels_per_m * quad->ref_width, quad->ipm_top);
	ipm_points[2] = Point2f(quad->ipm_left, quad->ipm_top +
			quad->pixels_per_m * quad->ref_height);
	ipm_points[3] = Point2f(quad->ipm_left +
			quad->pixels_per_m * quad->ref_width,
			quad->ipm_top +
			quad->pixels_per_m * quad->ref_height);

	ipm_matrix = getPerspectiveTransform(ref_points, ipm_points);

	hmat->a0 = ipm_matrix.at <double> (0, 0);
	hmat->a1 = ipm_matrix.at <double> (0, 1);
	hmat->a2 = ipm_matrix.at <double> (0, 2);
	hmat->b0 = ipm_matrix.at <double> (1, 0);
	hmat->b1 = ipm_matrix.at <double> (1, 1);
	hmat->b2 = ipm_matrix.at <double> (1, 2);
	hmat->c0 = ipm_matrix.at <double> (2, 0);
	hmat->c1 = ipm_matrix.at <double> (2, 1);
	hmat->c2 = ipm_matrix.at <double> (2, 2);

	return ipm_matrix;
}

static void apply_mat_to_image(const char *in, const char *out, Mat ipm_matrix)
{
	Mat src = imread(in);
	if (src.empty()) {
		printf("Could not open file %s\n", in);
		exit(1);
	}
	Mat dst = src.clone();

	warpPerspective(src, dst, ipm_matrix, dst.size());
	imwrite(out, dst);
}

static void dump_quadrangle(const struct quadrangle *cr)
{
	printf("\ttopl:\t%.1f, %.1f\n", cr->p_top_l.x, cr->p_top_l.y);
	printf("\ttopr:\t%.1f, %.1f\n", cr->p_top_r.x, cr->p_top_r.y);
	printf("\tbotl:\t%.1f, %.1f\n", cr->p_bot_l.x, cr->p_bot_l.y);
	printf("\tbotr:\t%.1f, %.1f\n", cr->p_bot_r.x, cr->p_bot_r.y);
	printf("\theight:\t%.1f\n", cr->ref_height);
	printf("\twidth:\t%.1f\n", cr->ref_width);
	printf("\tpixels:\t%.1f\n", cr->pixels_per_m);
	printf("\ttop:\t%.1f\n", cr->ipm_top);
	printf("\tleft:\t%.1f\n", cr->ipm_left);
}

static int read_rectangles_from_file(struct quadrangle *rect,
		const char *filename)
{
	FILE *file;
	struct quadrangle *r;

	file = fopen(filename, "r");
	if (!file) {
		fprintf(stderr, "File '%s' not found:", filename);
		perror("");
		return -1;
	}
	for (r = rect; r < (rect + 3); r++) {
		fscanf(file, "topl:\t%f, %f\n", &r->p_top_l.x, &r->p_top_l.y);
		fscanf(file, "topr:\t%f, %f\n", &r->p_top_r.x, &r->p_top_r.y);
		fscanf(file, "botl:\t%f, %f\n", &r->p_bot_l.x, &r->p_bot_l.y);
		fscanf(file, "botr:\t%f, %f\n", &r->p_bot_r.x, &r->p_bot_r.y);
		fscanf(file, "height:\t%f\n", &r->ref_height);
		fscanf(file, "width:\t%f\n", &r->ref_width);
		fscanf(file, "pixels:\t%f\n", &r->pixels_per_m);
		fscanf(file, "top:\t%f\n", &r->ipm_top);
		fscanf(file, "left:\t%f\n\n", &r->ipm_left);
		printf("\ndump rectangle 0x%X\n", r);
		dump_quadrangle(r);
	}

	fclose(file);
	return 0;
}

int main(int argc, const char** argv)
{
#define STR_SIZE 256
	static char rect_file[STR_SIZE];
	static char spl_file[STR_SIZE];
	static char rect1_file[STR_SIZE];
	static char rect2_file[STR_SIZE];
	static char rect3_file[STR_SIZE];
	static char matrix_file[STR_SIZE];
	struct hmatrix hmat;
	struct quadrangle rect[3];
	FILE *fhmatrix;
	Mat ipm_matrix;
	int ret;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <video path>\n", argv[0]);
		return 0;
	}

	snprintf(rect_file, STR_SIZE, "%s/olddata/rectangles.txt", argv[1]);
	snprintf(spl_file, STR_SIZE, "%s/olddata/sample.png", argv[1]);
	snprintf(rect1_file, STR_SIZE, "%s/olddata/rect1.png", argv[1]);
	snprintf(rect2_file, STR_SIZE, "%s/olddata/rect2.png", argv[1]);
	snprintf(rect3_file, STR_SIZE, "%s/olddata/rect3.png", argv[1]);
	snprintf(matrix_file, STR_SIZE, "%s/matrix.txt", argv[1]);

	ret = read_rectangles_from_file(rect, rect_file);
	if (ret) {
		return -1;
	}

	fhmatrix = fopen(matrix_file, "w");
	if (!fhmatrix) {
		fprintf(stderr, "Could not open file %s\n", rect_file);
		return -1;
	}

	ipm_matrix = generate_hmatrix(&hmat, &rect[0]);
	apply_mat_to_image(spl_file, rect1_file, ipm_matrix);
	print_hmatrix(&hmat, 1);
	save_hmatrix(fhmatrix, &hmat);
	printf("\n");

	ipm_matrix = generate_hmatrix(&hmat, &rect[1]);
	apply_mat_to_image(spl_file, rect2_file, ipm_matrix);
	print_hmatrix(&hmat, 2);
	save_hmatrix(fhmatrix, &hmat);
	printf("\n");

	ipm_matrix = generate_hmatrix(&hmat, &rect[2]);
	apply_mat_to_image(spl_file, rect3_file, ipm_matrix);
	print_hmatrix(&hmat, 3);
	save_hmatrix(fhmatrix, &hmat);
	printf("\n");

	fclose(fhmatrix);
	return 0;
}

