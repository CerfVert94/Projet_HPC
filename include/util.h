/* -------------- */
/* --- util.h --- */
/* -------------- */

#ifndef __UTIL_H__
#define __UTIL_H__

#define Vmin 1
#define Vmax 254

#define max(a,b) (a >= b ? a : b)
#define min(a,b) (a <= b ? a : b)

extern char *nom_func;

#define PASS(nom, test) fprintf(stderr, "%s%s    \t----- PASS -----\n", nom, test)
#define FAIL(nom, test) fprintf(stderr, "%s%s    \txxxxx FAIL xxxxx\n", nom, test); exit(EXIT_FAILURE)

#define error(msg) fprintf(stderr, "%s%s\n", msg, __func__); exit(EXIT_FAILURE)

#define UTEST(cond, test) if (cond) PASS(nom_func, test); else {FAIL(nom_func, test);}

#endif // __UTIL_H__