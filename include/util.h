/* -------------- */
/* --- util.h --- */
/* -------------- */

#ifndef __UTIL_H__
#define __UTIL_H__

#define Vmin 1
#define Vmax 254

#define BORD 2

#define max(a,b) (a >= b ? a : b)
#define min(a,b) (a <= b ? a : b)

// extern const char *nom_func;

#define PASS(func, test) fprintf(stderr, "%s%s   \t----- PASS -----\n", func, test)
#define FAIL(func, test) fprintf(stderr, "%s%s   \txxxxx FAIL xxxxx\n", func, test); exit(EXIT_FAILURE)

#define error(msg) fprintf(stderr, "%s%s\n", msg, __func__); exit(EXIT_FAILURE)

#define UNIT_TEST(cond, func, test) if (cond) PASS(func, test); else {FAIL(func, test);}

#define TIME(v)   


#define LALIGNED_STR     "%-45s"
#define RALIGNED_STR     "%45s"
#define RALIGNED_SINT    "%40d"
#define RALIGNED_UINT    "%40u"
#define RALIGNED_SLONG   "%40ld"
#define RALIGNED_ULONG   "%40lu"
#define RALIGNED_SLLONG  "%40lld"
#define RALIGNED_ULLONG  "%40llu"
#define RALIGNED_FLOAT   "%40.2f"
#define RALIGNED_DOUBLE  "%40.2lf"


static inline void exit_on_error(const char *msg){
    perror(msg);
    exit(EXIT_FAILURE);
}
static void exit_on_error(const char *msg);

void binary_to_octal_ui8matrix(uint8 **ppInput, long nrl, long nrh, long ncl, long nch);
void octal_to_binary_ui8matrix(uint8 **ppInput, long nrl, long nrh, long ncl, long nch);


#define NROW(nrl, nrh)  (nrh - nrl + 1)
#define NCOL(ncl, nch)  (nch - ncl + 1)
#endif // __UTIL_H__