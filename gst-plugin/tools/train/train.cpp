#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "xml/pugixml.hpp"
#include "svm/svm.h"
#include "thog/thog.h"
#include "png.h"

static void usage(char *progname)
{
	printf("Usage: %s\n"
		"\t\t<input file> e.g. icdar/positives.txt\n"
		"\t\t<output file> e.g. out/positives.txt\n"
		"\t\t<label> e.g. 1\n"
		"\t\t<thog settings> e.g. input/1_7_9.txt\n", progname);
}

int main(int argc, char *argv[])
{
	FILE *fin;
	FILE *fout;
	int label;
	char *thog_settings;

	/* Check input parameters: */
	if (argc != 5) {
		usage(argv[0]);
		return -1;
	}

	fin = fopen(argv[1], "r");
	if (NULL == fin) {
		perror("train");
		return -1;
	}
	fout = fopen(argv[2], "w"); 
	if (NULL == fout) {
		perror("train");
		fclose(fin);
		return -1;
	}
	label = atoi(argv[3]); 
	thog_settings = argv[4];

	/* Loading T-HOG settings: */
	struct_thog sthog = load_settings(thog_settings);

	while (!feof(fin)) {
		char *filename = (char *) malloc(512 * sizeof(char));
		unsigned error;
		unsigned char* colorpng;
		unsigned int nrows;
		unsigned int ncols;
		int match;

		match = fscanf(fin, "%s", filename);
		if (!feof(fin) && match) {
			error = lodepng_decode32_file(&colorpng, &ncols,
					&nrows, filename);
			if (error) {
				printf("error %u: %s\n", error,
						lodepng_error_text(error));
			}
			unsigned char* image =
				convert_argb_to_gray(colorpng, nrows, ncols);
			free(colorpng);
			double *descriptor = thog(image, nrows, ncols, sthog);
			free(image);
			int i;
			fprintf(fout, "%d ", label);
			for (i = 0; i < sthog.nob; i++) {
				fprintf(fout, "%d:%f ", i + 1, descriptor[i]);
			}
			fprintf(fout, "\n");
		}
		free(filename);
	}
	fclose(fin);
	fclose(fout);

	return 0;
}

