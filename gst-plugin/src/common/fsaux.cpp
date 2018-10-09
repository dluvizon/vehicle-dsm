#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "fsaux.h"
#include "debug.h"

/* 
 * thaks for
 * http://nion.modprobe.de/blog/archives/357-Recursive-directory-creation.html
 */
int fs_mkdir_recursive(const char *dir)
{
	char tmp[256];
	char *p = NULL;
	size_t len;
	int ret;

	assert(dir);

	ret = 0;
	len = strlen(dir);
	if (len >= sizeof(tmp)) {
		print_err("in _%s_ directory path too large", __func__);
		ret = -1;
		goto end;
	}
	snprintf(tmp, sizeof(tmp), "%s", dir);
	if(tmp[len - 1] == '/') {
		/* remove trailing slash */
		tmp[len - 1] = '\0';
	}

	for(p = tmp + 1; *p; p++) {
		if (*p == '/') {
			*p = '\0';
			if (!fs_dir_exist(tmp)) {
				ret = mkdir(tmp, S_IRWXU);
				if (-1 == ret) {
					goto end;
				}
			}
			*p = '/';
		}
	}

	if (!fs_dir_exist(tmp)) {
		/* if this dir is not here, create it */
		ret = mkdir(tmp, S_IRWXU);
	}

end:
	return ret;
}

bool fs_dir_exist(const char *dir)
{
	bool ret;
	DIR *ds;

	assert(dir);

	ds = opendir(dir);
	ret = !!ds;
	if (ret) {
		/* if directory exist, close it */
		closedir(ds);
	}

	return ret;
}
