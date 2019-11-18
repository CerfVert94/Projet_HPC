/* -------------- */
/* --- util.h --- */
/* -------------- */

#ifndef __UTIL_H__
#define __UTIL_H__

#define Vmin 1
#define Vmax 254

#define max(a,b) (a >= b ? a : b)
#define min(a,b) (a <= b ? a : b)

// extern const char *nom_func;

#define PASS(func, test) fprintf(stderr, "%s%s   \t----- PASS -----\n", func, test)
#define FAIL(func, test) fprintf(stderr, "%s%s   \txxxxx FAIL xxxxx\n", func, test); exit(EXIT_FAILURE)

#define error(msg) fprintf(stderr, "%s%s\n", msg, __func__); exit(EXIT_FAILURE)

#define UNIT_TEST(cond, func, test) if (cond) PASS(func, test); else {FAIL(func, test);}

#endif // __UTIL_H__