#ifndef __STRAUX_H_
#define __STRAUX_H_

/**
 * Format the input string with %05d of val formated.
 * @param txt Input string.
 * @param val Value.
 * @return Formated string.
 */
char *format_uint_in_5x(char *txt, unsigned int val);

/**
 * Replace the string pat by rep in str string.
 */
char *replace_string(char *str, const char *pat, const char *rep);

#endif /* __STRAUX_H_  */
