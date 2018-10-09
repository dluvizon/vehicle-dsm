/**
 * error.h
 */

#ifndef __ERROR_H
#define __ERROR_H

#include <stdio.h>
#include <stdarg.h>

void KLTError(const char *fmt, ...);
void KLTWarning(const char *fmt, ...);

#endif /* __ERROR_H */

