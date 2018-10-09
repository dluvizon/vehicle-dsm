#ifndef __FONT_H
#define __FONT_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void draw_text(unsigned char *image, int width, int x, int y, int scale,
		int font_color, int back_color, const char *msg, ...);

#endif /* __FONT_H */
