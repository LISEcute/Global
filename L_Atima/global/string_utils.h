#ifndef __STRING_UTILS
#define __STRING_UTILS

#include <ctype.h>
/*
 * Tong Zhang <zhangt@frib.msu.edu>
 * 2020-08-27 09:24:28 EDT
 */

char *strlwrLoc(char *s)
{
    char *str_lower = s;
    int i = 0;
    while (s[i]) {
        str_lower[i] = tolower(s[i]);
        i++;
    }
    return str_lower;
}

char *struprLoc(char *s)
{
    char *str_upper = s;
    int i = 0;
    while (s[i]) {
        str_upper[i] = toupper(s[i]);
        i++;
    }
    return str_upper;
}

#endif // __STRING_UTILS
