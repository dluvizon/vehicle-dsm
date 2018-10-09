#ifndef __FSAUX_H_
#define __FSAUX_H_

#include <stdbool.h>

/**
 * Create a new directory recursively, if needed.
 * @param dir Path to directory.
 */
int fs_mkdir_recursive(const char *dir);

/**
 * Return true if this directory exist, false if not.
 * @param dir Path to directory.
 * @return True or false.
 */
bool fs_dir_exist(const char *dir);

#endif /* __FSAUX_H_  */
