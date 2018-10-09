/**
 * @file swt-wrapper.h
 * @author Diogo Luvizon <diogo at luvizon dot com>
 * @date 28/04/2015
 */
#ifndef __SWT_WRAPPER_H
#define __SWT_WRAPPER_H

#include <vector>

#include "backgroundsub/backgroundsub.h"
#include "snoopertext/snoopertext.h"
#include "seg/toggle.h"
#include "svm/svm.h"
#include "thog/thog.h"
#include "image.h"

vector<plate> detection_swt(
                unsigned char *image,
		struct snoopertext *s,
                int ncols, int nrows,
                char *out_path,
                vector<struct sub_slope> &slopes,
                int margin_X, int margin_Y,
                int iframe);

#endif /* __SWT_WRAPPER_H */
