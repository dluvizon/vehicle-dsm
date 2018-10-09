/**
 * @file detection-minetto.h
 * @author Rodrigo Minetto
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 10/03/2015
 */

#ifndef __DETECTION_MINETTO_H
#define __DETECTION_MINETTO_H

#include <vector>

#include "backgroundsub/backgroundsub.h"
#include "seg/toggle.h"
#include "svm/svm.h"
#include "thog/thog.h"
#include "snoopertext.h"
#include "image.h"

#define DETECT_MOTOCYCLE 1
#define CAR_MIN_WIDTH 320

vector<plate> detection_minetto(
                unsigned char *image,
                struct snoopertext *s,
                int ncols, int nrows,
                char *out_path,
                char *out_image_name,
                vector<struct sub_slope> &slopes,
                int margin_X, int margin_Y,
                int iframe,
                unsigned char *umask);

inline bool slope_is_moto(struct sub_slope s)
{
	/* If not detecting motocycles, just return false. */
#if !DETECT_MOTOCYCLE
	return false;
#else
	/* Check slope consistency. */
	if (s.right < s.left) {
		return false;
	}
	/* If is larger than it, then is a car. */
	if (s.right - s.left > CAR_MIN_WIDTH) {
		return false;
	}
	/* Otherwise, it must be a motocycle. */
	return true;
#endif
}

#endif /* __DETECTION_MINETTO_H */
