/* -------------- */
/* --- util.h --- */
/* -------------- */

#ifndef __UTIL_H__
#define __UTIL_H__

#define Vmin 1
#define Vmax 254

#define max(a,b) (a >= b ? a : b)
#define min(a,b) (a <= b ? a : b)

#define PASS() fprintf(stderr, "%s\t=\t----- PASS -----\n", __func__)
#define FAIL() fprintf(stderr, "%s\t=\txxxxx FAIL xxxxx\n", __func__); exit(EXIT_FAILURE)

#endif // __UTIL_H__