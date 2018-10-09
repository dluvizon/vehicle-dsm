/**
 * @file detection-zheng.h
 * @author Rodrigo Minetto
 * @date 05/05/2015
 */

#ifndef __DETECTION_ZHENG_H
#define __DETECTION_ZHENG_H

#include <vector>

#include "backgroundsub/backgroundsub.h"
#include "seg/toggle.h"
#include "svm/svm.h"
#include "thog/thog.h"
#include "snoopertext.h"
#include "image.h"

vector<plate> detection_zheng(
                unsigned char *image,
                struct snoopertext *s,
                int ncols, int nrows,
                char *out_path,
                char *out_image_name,
                vector<struct sub_slope> &slopes,
                int margin_X, int margin_Y,
                int iframe,
                unsigned char *umask);

#endif /* __DETECTION_ZHENG_H */
