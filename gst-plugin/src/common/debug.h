/** LightWeight Computer Vision Library
 * ------------------------------------
 * @Copyright (C) 2014 SVi Smart Vision
 * @file debug.h
 * @brief Defines some macros and inline functions for debug pourpose.
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 23/10/2014
 *
 * It is strongly recommended to not include this file in *.h files.
 * Use this from *.c or *.cpp files and define the macro _DLEVEL as desired.
 */

#ifndef __LWCV_DEBUG_H
#define __LWCV_DEBUG_H

/* set default debug state if not defined */
#ifndef _DLEVEL
#define _DLEVEL 1
#endif

#include <sys/timeb.h>
#include <stdarg.h>
#include <stdio.h>
#include <time.h>

#include "term-colors.h"
#include "config.h"

#define _PACK_INFO_STR \
	"Package: " cterm(TB_WHT, PACKAGE_NAME) " version " PACKAGE_VERSION \
	" in " PACKAGE_DATE "\nBug Report: " PACKAGE_BUGREPORT

#undef USE_TIMESTAMP

/**
 * The variable 'volatile int _dbg_tab' must be defined in just ONE location,
 * inside a *.c file.  This variable will be used by print_tab function.
 */
extern volatile int _dbg_tab;
const int _tab_max = 10;

enum _ptab {
	DBG_CALL,
	DBG_EXIT,
};

inline void print_err(const char *msg, ...)
{
	va_list argptr;

	va_start(argptr, msg);
	fprintf(stderr, cterm(TB_RED, "error") ": ");
	vfprintf(stderr, msg, argptr);
	fprintf(stderr, "\n" _PACK_INFO_STR "\n\n");
	va_end(argptr);
}

#if _DLEVEL > 1
inline void print_tab(enum _ptab t, const char *msg, ...)
{
	va_list argptr;
	int i;

	va_start(argptr, msg);
	fprintf(stdout, cterm(TB_WHT, PACKAGE_NAME) ":");
	for (i = 0; (i < _dbg_tab) && (i < _tab_max); i++) {
		fprintf(stdout, "  ");
	}
	if (DBG_CALL == t) {
		fprintf(stdout, "  » ");
		_dbg_tab++;
	} else if (DBG_EXIT == t) {
		fprintf(stdout, "« ");
		_dbg_tab--;
	} else {
		fprintf(stdout, "? ");
	}
	if (_dbg_tab > _tab_max) {
		_dbg_tab = _tab_max;
	}
	if (_dbg_tab < 0) {
		_dbg_tab = 0;
	}
	vfprintf(stdout, msg, argptr);
	va_end(argptr);
}
#else
# define print_tab(...)
#endif


#if _DLEVEL > 1
inline void print_dbg(const char *msg, ...)
{
	va_list argptr;

	va_start(argptr, msg);
	fprintf(stdout, cterm(TB_WHT, PACKAGE_NAME) ": "
			cterm(TR_BLU, "debug") ": ");
	vfprintf(stdout, msg, argptr);
	va_end(argptr);
}
#else
# define print_dbg(...)
#endif

#if _DLEVEL > 0
inline void print_warn(const char *msg, ...)
{
	va_list argptr;

	va_start(argptr, msg);
	fprintf(stderr, cterm(TB_WHT, PACKAGE_NAME) ": "
			cterm(TB_YLW, "warning") ": ");
	vfprintf(stderr, msg, argptr);
	va_end(argptr);
}
#else
# define print_warn(...)
#endif

#if _DLEVEL > 0
inline void print_timestamp(const char *msg, ...)
{
	va_list argptr;

#ifdef USE_TIMESTAMP
	struct timeb tb;

	ftime(&tb);
	va_start(argptr, msg);
	fprintf(stdout, cterm(TB_WHT, PACKAGE_NAME) ": " cterm(TB_GRN, "ts")
			": %.8s-%03d: ", ctime(&tb.time) + 11, tb.millitm);
#else
	va_start(argptr, msg);
	fprintf(stdout, cterm(TB_WHT, PACKAGE_NAME) ": " cterm(TB_GRN, "ts")
			": 00:00:00-000: ");
#endif
	vfprintf(stdout, msg, argptr);
	va_end(argptr);
}
#else
# define print_timestamp(...)
#endif

#endif /* __LWCV_DEBUG_H */
