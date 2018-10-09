#ifndef __SPEED_H
#define __SPEED_H

#include <stdio.h>

#include "hmatrix.h"
#include "klt/klt.h"

#define DEBUG_OUTLIER_KLT 0

typedef struct _speed_i {
	double v_x;
	double v_y;
	double m_x;
	double m_y;
} speed_i;

int load_ipm_matrix_from_file(const char *filename);

void compute_ipm_features(KLT_FeatureList features, int iframe);

speed_i compute_velocity_vector(KLT_FeatureList features,
		int frame_i, int frame_j);

double metric_calc_velocity_ipm(speed_i speed, unsigned lane);

void outlier_removal(uint8_t *image, int width, int height,
		KLT_FeatureList features, int iframe);

#endif /* __SPEED_H */
