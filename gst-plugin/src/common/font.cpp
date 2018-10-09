/**
 * @file cricfont.c
 * @brief A simple 7-bit ASCII font for use with EricDraw
 * @author Eric Shalov Jan 2010
 */

#include <stdarg.h>

#include "font.h"

const static char font[128][72 + 1] = {
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X X X X "
	"X     X "
	"X XXX X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	" XXXXX  "
	"XXXXXXX "
	"XX X XX "
	"XXXXXXX "
	"XX   XX "
	"XXXXXXX "
	" XXXXX  "
	"        ", 

	"        "
	"        "
	" XX XX  "
	"XXXXXXX "
	"XXXXXXX "
	" XXXXX  "
	"  XXX   "
	"   X    "
	"        ", 

	"        "
	"   X    "
	"  XXX   "
	" XXXXX  "
	"XXXXXXX "
	" XXXXX  "
	"  XXX   "
	"   X    "
	"        ", 

	"        "
	"  XXX   "
	"X  X  X "
	"XXXXXXX "
	"X  X  X "
	"  XXX   "
	"   X    "
	"  XXX   "
	"        ", 

	"        "
	"   X    "
	"  XXX   "
	" XXXXX  "
	"XXXXXXX "
	" XXXXX  "
	"   X    "
	" XXXXX  "
	"        ", 

	"        "
	"        "
	"        "
	"  XXX   "
	" XXXXX  "
	"  XXX   "
	"        "
	"        "
	"        ", 

	"        "
	"XXXXXXX "
	"XXXXXXX "
	"XX   XX "
	"X     X "
	"XX   XX "
	"XXXXXXX "
	"XXXXXXX "
	"        ", 

	"        "
	"        "
	"        "
	"  XXX   "
	" XX XX  "
	"  XXX   "
	"        "
	"        "
	"        ", 

	"        "
	"XXXXXXX "
	"XXXXXXX "
	"XX   XX "
	"X  X  X "
	"XX   XX "
	"XXXXXXX "
	"XXXXXXX "
	"        ", 

	"        "
	"    XXX "
	"     XX "
	"    X X "
	"   X    "
	"  XXX   "
	" X   X  "
	"  XXX   "
	"        ", 

	"        "
	"  XXX   "
	" X   X  "
	"  XXX   "
	"   X    "
	" XXXXX  "
	"   X    "
	"   X    "
	"        ", 

	"        "
	"  XXXX  "
	"  X  X  "
	"  XXXX  "
	"  X     "
	" XX     "
	"XXX     "
	" XX     "
	"        ", 

	"        "
	"  XXXXX "
	"  X   X "
	"  XXXXX "
	"  X   X "
	" XX  XX "
	"XXX XXX "
	" XX  XX "
	"        ", 

	"        "
	" X X X  "
	"  XXX   "
	" X X X  "
	" XXXXX  "
	" X X X  "
	"  XXX   "
	" X X X  "
	"        ", 

	"        "
	"  XX    "
	"  XXX   "
	"  XXXX  "
	"  XXXXX "
	"  XXXX  "
	"  XXX   "
	"  XX    "
	"        ", 

	"        "
	"     XX "
	"    XXX "
	"   XXXX "
	"  XXXXX "
	"   XXXX "
	"    XXX "
	"     XX "
	"        ", 

	"        "
	"  XX    "
	" XXXX   "
	"XXXXXX  "
	"  XX    "
	"XXXXXX  "
	" XXXX   "
	"  XX    "
	"        ", 

	"        "
	" X X    "
	" X X    "
	" X X    "
	" X X    "
	" X X    "
	"        "
	" X X    "
	"        ", 

	"        "
	"  XXXXX "
	" X  X X "
	" X  X X "
	"  XXX X "
	"    X X "
	"    X X "
	"    X X "
	"        ", 

	"        "
	"  XXX   "
	" XX     "
	"  XXX   "
	"  X X   "
	"  XXX   "
	"     XX "
	"   XXX  "
	"        ", 

	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"XXXXXXX "
	"XXXXXXX "
	"        ", 

	"        "
	"   X    "
	"  XXX   "
	" X X X  "
	"   X    "
	" X X X  "
	"  XXX   "
	"   X    "
	"XXXXXXX ", 

	"        "
	"   X    "
	"  XXX   "
	" X X X  "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"        ", 

	"        "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	" X X X  "
	"  XXX   "
	"   X    "
	"        ", 

	"        "
	"        "
	"  XX    "
	"   XX   "
	"XXXXXX  "
	"   XX   "
	"  XX    "
	"        "
	"        ", 

	"        "
	"        "
	"  XX    "
	" XX     "
	"XXXXXX  "
	" XX     "
	"  XX    "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"X       "
	"X       "
	"XXXXXX  "
	"        "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	" X   X  "
	"XXXXXXX "
	" X   X  "
	"        "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"   X    "
	"  XXX   "
	" XXXXX  "
	"XXXXXXX "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"XXXXXXX "
	" XXXXX  "
	"  XXX   "
	"   X    "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        ", 

	"        "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"        "
	"   X    "
	"        ", 

	"        "
	"  X X   "
	"  X X   "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        ", 

	"        "
	"  X  X  "
	" XXXXXX "
	"  X  X  "
	"  X  X  "
	" XXXXXX "
	"  X  X  "
	"        "
	"        ", 

	"        "
	" XXXXX  "
	"X  X    "
	"X  X    "
	" XXXXX  "
	"   X  X "
	"   X  X "
	" XXXXX  "
	"   X    ", 

	"        "
	"  XX    "
	" X  X X "
	"  XX X  "
	"    X   "
	"   X    "
	"  X XX  "
	" X X  X "
	"X   XX  ", 

	"        "
	" XXX    "
	"X   X   "
	" X X    "
	"  X     "
	" X X X  "
	"X   X   "
	" XXX X  "
	"        ", 

	"        "
	"   X    "
	"   X    "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        ", 

	"        "
	"   X    "
	" XX     "
	"X       "
	"X       "
	"X       "
	" XX     "
	"   X    "
	"        ", 

	"        "
	"   X    "
	"    XX  "
	"      X "
	"      X "
	"      X "
	"    XX  "
	"   X    "
	"        ", 

	"        "
	"   X    "
	" X X X  "
	"  XXX   "
	"   X    "
	"  XXX   "
	" X X X  "
	"   X    "
	"        ", 

	"        "
	"        "
	"   X    "
	"   X    "
	" XXXXX  "
	"   X    "
	"   X    "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"   X    "
	"   XX   "
	"    X   ", 

	"        "
	"        "
	"        "
	"        "
	"XXXXXXX "
	"        "
	"        "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"X       "
	"        ", 

	"        "
	"      X "
	"     X  "
	"    X   "
	"   X    "
	"  X     "
	" X      "
	"X       "
	"        ", 

	"        "
	"  XXX   "
	" X   X  "
	"X   X X "
	"X  X  X "
	"X X   X "
	" X   X  "
	"  XXX   "
	"        ", 

	"        "
	"   X    "
	"  XX    "
	" X X    "
	"   X    "
	"   X    "
	"   X    "
	" XXXXX  "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"     XX "
	"   XX   "
	"  X     "
	" X      "
	"XXXXXXX "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"      X "
	"   XXX  "
	"      X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	"   XXXX "
	"  X   X "
	" X    X "
	"XXXXXXX "
	"      X "
	"      X "
	"      X "
	"        ", 

	"        "
	"XXXXXXX "
	"X       "
	"X       "
	"XXXXXX  "
	"      X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	" XXXXX  "
	"X       "
	"X       "
	"XXXXXX  "
	"X     X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	"XXXXXXX "
	"      X "
	"     X  "
	"    X   "
	"   X    "
	"  X     "
	" X      "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X     X "
	" XXXXX  "
	"X     X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X     X "
	" XXXXX  "
	"      X "
	"      X "
	" XXXXX  "
	"        ", 

	"        "
	"        "
	"  X     "
	"        "
	"        "
	"        "
	"  X     "
	"        "
	"        ", 

	"        "
	"        "
	"  X     "
	"        "
	"        "
	"        "
	"  XX    "
	"   X    "
	"        ", 

	"        "
	"    XX  "
	"   XX   "
	"  XX    "
	" XX     "
	"  XX    "
	"   XX   "
	"    XX  "
	"        ", 

	"        "
	"        "
	"        "
	"XXXXXXX "
	"        "
	"XXXXXXX "
	"        "
	"        "
	"        ", 

	"        "
	"XX      "
	" XX     "
	"  XX    "
	"   XX   "
	"  XX    "
	" XX     "
	"XX      "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X    XX "
	"    XX  "
	"        "
	"   XX   "
	"   XX   "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X XXX X "
	"X X X X "
	"X XXXXX "
	"X       "
	" XXXXX  "
	"        ", 

	"        "
	"   X    "
	"  X X   "
	" X   X  "
	"X     X "
	"XXXXXXX "
	"X     X "
	"X     X "
	"        ", 

	"        "
	"XXXXXX  "
	"X     X "
	"X     X "
	"XXXXXX  "
	"X     X "
	"X     X "
	"XXXXXX  "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X       "
	"X       "
	"X       "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	"XXXXXX  "
	"X     X "
	"X     X "
	"X     X "
	"X     X "
	"X     X "
	"XXXXXX  "
	"        ", 

	"        "
	"XXXXXXX "
	"X       "
	"X       "
	"XXXXXXX "
	"X       "
	"X       "
	"XXXXXXX "
	"        ", 

	"        "
	"XXXXXXX "
	"X       "
	"X       "
	"XXXXXXX "
	"X       "
	"X       "
	"X       "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X       "
	"X  XXX  "
	"X     X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	"X     X "
	"X     X "
	"X     X "
	"XXXXXXX "
	"X     X "
	"X     X "
	"X     X "
	"        ", 

	"        "
	"XXXXXXX "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"XXXXXXX "
	"        ", 

	"        "
	"XXXXXXX "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"X  X    "
	" XX     "
	"        ", 

	"        "
	"X    X  "
	"X   X   "
	"X  X    "
	"X XX    "
	"XX  X   "
	"X    X  "
	"X     X "
	"        ", 

	"        "
	"X       "
	"X       "
	"X       "
	"X       "
	"X       "
	"X       "
	"XXXXXXX "
	"        ", 

	"        "
	"XX   XX "
	"X X X X "
	"X  X  X "
	"X  X  X "
	"X     X "
	"X     X "
	"X     X "
	"        ", 

	"        "
	"X     X "
	"XXX   X "
	"X XX  X "
	"X  XX X "
	"X   XXX "
	"X    XX "
	"X     X "
	"        ", 

	"        "
	"  XXX   "
	" X   X  "
	"X     X "
	"X     X "
	"X     X "
	" X   X  "
	"  XXX   "
	"        ", 

	"        "
	"XXXXXX  "
	"X     X "
	"X     X "
	"XXXXXX  "
	"X       "
	"X       "
	"X       "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X     X "
	"X     X "
	"X   X X "
	"X    X  "
	" XXXX X "
	"        ", 

	"        "
	"XXXXXX  "
	"X     X "
	"X     X "
	"XXXXXX  "
	"X     X "
	"X     X "
	"X     X "
	"        ", 

	"        "
	" XXXXX  "
	"X     X "
	"X       "
	" XXXXX  "
	"      X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	"XXXXXXX "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"   X    "
	"        ", 

	"        "
	"X     X "
	"X     X "
	"X     X "
	"X     X "
	"X     X "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	"X     X "
	"X     X "
	"X     X "
	"X     X "
	" X   X  "
	"  X X   "
	"   X    "
	"        ", 

	"        "
	"X     X "
	"X     X "
	"X  X  X "
	"X  X  X "
	"X  X  X "
	"X  X  X "
	" XX XX  "
	"        ", 

	"        "
	"X     X "
	" X   X  "
	"  X X   "
	"   X    "
	"  X X   "
	" X   X  "
	"X     X "
	"        ", 

	"        "
	"X     X "
	"X     X "
	" X   X  "
	"  X X   "
	"   X    "
	"   X    "
	"   X    "
	"        ", 

	"        "
	"XXXXXXX "
	"     X  "
	"    X   "
	"   X    "
	"  X     "
	" X      "
	"XXXXXXX "
	"        ", 

	"        "
	" XXXXX  "
	" X      "
	" X      "
	" X      "
	" X      "
	" X      "
	" XXXXX  "
	"        ", 

	"        "
	"X       "
	" X      "
	"  X     "
	"   X    "
	"    X   "
	"     X  "
	"      X "
	"        ",

	"        "
	" XXXXX  "
	"     X  "
	"     X  "
	"     X  "
	"     X  "
	"     X  "
	" XXXXX  "
	"        ", 

	"        "
	"   X    "
	"  X X   "
	" X   X  "
	"        "
	"        "
	"        "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        "
	"XXXXXXXX", 

	"        "
	"  XX    "
	"   XX   "
	"        "
	"        "
	"        "
	"        "
	"        "
	"        ", 

	"        "
	"        "
	"        "
	"  XXX   "
	" X   X  "
	"X     X " 
	" X   XX "
	"  XXX X " 
	"        ", 

	"        "
	"X       "
	"X       "  
	"X XXX   "
	"XX   X  "
	"X     X "
	"X     X "  
	"XXXXXX  "
	"        ",

	"        "
	"        "
	"        "
	" XXXXX  "
	"X     X "
	"X       "
	"X     X "
	" XXXXX  "
	"        ", 

	"        "
	"      X "
	"      X "
	"  XXX X "
	"XX   XX "
	"X     X "
	"X     X "
	" XXXXXX "
	"        ", 

	"        "
	"        "
	"        "
	" XXXXX  "
	"X     X "
	"XXXXXX  "
	"X       "
	" XXXXX  "
	"        ", 

	"        "
	"  XXXX  "
	" X    X "
	" X      "
	"XXXXXX  "
	" X      "
	" X      "
	" X      "
	"        ", 

	"        "
	"        "
	"        "
	"  XXXX  "
	" X    X "
	"  XXXXX "
	"      X "
	" X   X  "
	"  XXX   ", 

	"        "
	"X       "
	"X       "
	"X       "
	"X XXX   "
	"XX   X  "
	"X     X "
	"X     X "
	"        ", 

	"        "
	"        "
	"  X     "
	"        "
	"  X     "
	"  X     "
	"  X     "
	"  X     "
	"        ", 

	"        "
	"        "
	"   X    "
	"        "
	"   X    "
	"   X    "
	"X  X    "
	"X  X    "
	" XX     ", 

	"        "
	"        "
	"X       "
	"X  X    "
	"X X     "
	"XXX     "
	"X  X    "
	"X   X   "
	"        ", 

	"        "
	" X      "
	" X      "
	" X      "
	" X      "
	" X      "
	" X      "
	" XX     "
	"        ", 

	"        "
	"        "
	"        "
	"X XXXX  "
	"XX X XX "
	"X  X  X "
	"X  X  X "
	"X  X  X " 
	"        ", 

	"        "
	"        "
	"        "
	"X XXX   "
	"XX  XX  "
	"X    X  "
	"X    X  "
	"X    X  "
	"        ", 

	"        "
	"        "
	"        "
	" XXXX   "
	"X    X  "
	"X    X  "
	"X    X  "
	" XXXX   "
	"        ", 

	"        "
	"        "
	"X       "
	"X XX    "
	"XX  X   "
	"X   X   "
	"XXXX    "
	"X       "
	"X       ", 

	"        "
	"        "
	"        "
	"  XXXX  "
	" X   X  "
	"X    X  "
	"X    X  "
	" XXXXX  "
	"     X  ", 

	"        "
	"        "
	"X       "
	"XXXXX   "
	"X    X  "
	"X       "
	"X       "
	"X       "
	"        ", 

	"        "
	"        "
	"        "
	" XXXXX  "
	"X       "
	" XXXXX  "
	"      X "
	" XXXXX  "
	"        ", 

	"        "
	" X      "
	" X      "
	" X      "
	"XXXX    "
	" X      "
	" X   X  "
	"  XXX   "
	"        ", 

	"        "
	"        "
	"        "
	"X    X  "
	"X    X  "
	"X    X  "
	"X   XX  "
	" XXX X  "
	"        ", 

	"        "
	"        "
	"        "
	"X   X   "
	"X   X   "
	"X   X   "
	" X X    "
	"  X     "
	"        ", 

	"        "
	"        "
	"        "
	"X     X "
	"X  X  X "
	"X  X  X "
	"X  X  X "
	" XX XX  "
	"        ", 

	"        "
	"        "
	"        "
	"X   X   "
	" X X    "
	"  X     "
	" X X    "
	"X   X   "
	"        ", 

	"        "
	"        "
	"        "
	" X   X  "
	" X   X  "
	"  XXXX  "
	"     X  "
	" X   X  "
	"  XXX   ", 

	"        "
	"        "
	"        "
	"XXXXX   "
	"   X    "
	"  X     "
	" X      "
	"XXXXX   "
	"        ",

	"        "
	"  XX    "
	" X      "
	" X      "
	"X       "
	" X      "
	" X      "
	"  XX    "
	"        ", 

	"        "
	"  XX    "
	"  XX    "
	"  XX    "
	"  XX    "
	"  XX    "
	"  XX    "
	"  XX    "
	"        ", 

	"        "
	" XX     "
	"   X    "
	"   X    "
	"    X   "
	"   X    "
	"   X    "
	" XX     "
	"        ", 

	"        "
	"        "
	"        "
	" XX     "
	"X  X  X "
	"    XX  "
	"        "
	"        "
	"        ", 

	"        "
	"   X    "
	"  X X   "
	" X   X  "
	"X     X "
	"X     X "
	"XXXXXXX "
	"        "
	"        "
};

void draw_text(unsigned char *image, int width, int x, int y, int scale,
		int font_color, int back_color, const char *msg, ...)
{
	va_list argptr;
	char text[256];
	int p;
	int n;
	int h;
	int v;

	va_start(argptr, msg);
	va_end(argptr);
	vsnprintf(text, sizeof(text), msg, argptr);
	va_end(argptr);

	for (n = 0; text[n]; n++) {
		for (h = 0; h < 8; h++) {
			for (v = 0; v < 9; v++) {
				if (scale == 1) {
					p = (y + v) * width + (x + h + n * 8);
					image[p] = font[(unsigned char)
						text[n]][v * 9 + h] == 'X' ?
						255 : 0;
				} else {
					int x1 = x + n * 8 * scale + h * scale;
					int y1 = y + v * scale;
					int x2 = x + n * 8 * scale +
						h * scale + (scale - 1);
					int y2 = y + v * scale + (scale - 1);
					for (int yy = y1; yy <= y2; yy++) {
						for (int xx = x1; xx <=
								x2; xx++) {
							p = yy * width + xx;
							image[p] = (font[(unsigned char) text[n]][v * 8 + h] & 0x7F) == 'X' ? font_color : back_color;
						}
					}
				}
			}
		}
	}
}

