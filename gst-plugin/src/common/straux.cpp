#include <string.h>
#include <stdio.h>

#include "straux.h"

char *format_uint_in_5x(char *txt, unsigned int val)
{
	char strval[8];

	sprintf(strval, "%05d", val);

	return replace_string(txt, "xxxxx", strval);
}

char *replace_string(char *str, const char *pat, const char *rep)
{
	int ret;
	char *p;

	/* is 'pat' even in 'str'? */
	p = strstr(str, pat);
	if (!p) {
		return NULL;
	}
	strncpy(p, rep, strlen(pat));

	return str;
}
