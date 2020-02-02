/* -------------- */
/* --- util.h --- */
/* -------------- */

#ifndef __UTIL_H__
#define __UTIL_H__

typedef enum {SCALAR, SIMD} instruction_type;
#define Vmin 1
#define Vmax 254

#define BORD 2
#define vBORD (BORD/16 +1)

#define max(a,b) (a >= b ? a : b)
#define min(a,b) (a <= b ? a : b)

// extern const char *nom_func;

#define PASS(func, test) fprintf(stderr, "%s%s   \t----- PASS -----\n", func, test)
#define FAIL(func, test) fprintf(stderr, "%s%s   \txxxxx FAIL xxxxx\n", func, test); exit(EXIT_FAILURE)

#define error(msg) fprintf(stderr, "%s%s\n", msg, __func__); exit(EXIT_FAILURE)

#define UNIT_TEST(cond, func, test) if (cond) PASS(func, test); else {FAIL(func, test);}

#define TIME(v)   


#define LALIGNED_STR     "%-60s"
#define RALIGNED_STR     "%60s"
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


typedef uint8** (*pack_func_t)(uint8 **X, long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch);
typedef uint8** (*unpack_func_t)(uint8 **X, long nrl, long nrh, long ncl, long nch);

uint8 **fcpacked_ui8matrix (           long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord);
uint8 **frpacked_ui8matrix (           long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord);
uint8 **hcpacked_ui8matrix (           long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord);
uint8 **hrpacked_ui8matrix (           long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord);
void free_packed_ui8matrix (uint8 **X                                        , long  packed_nrl, long  packed_nrh, long  packed_ncl, long  packed_nch, long  bord);

void   fcpack_ui8matrix_ui8matrix (uint8 **X, long nrl, long nrh, long ncl, long nch, long packed_nrl, long packed_nrh, long packed_ncl, long packed_nch, long bord, uint8 **Y);
void unfcpack_ui8matrix_ui8matrix (uint8 **X, long nrl, long nrh, long ncl, long nch, long packed_nrl, long packed_nrh, long packed_ncl, long packed_nch, long bord, uint8 **Y);
// uint8 **unhpack_binary_ui8matrix(uint8 **X, long nrl, long nrh, long ncl, long nch);

uint8 **vpack_binary_ui8matrix  (uint8 **X, long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch);
uint8 **unvpack_binary_ui8matrix(uint8 **X, long nrl, long nrh, long ncl, long nch);

static inline long roundup_over8(long n);
static inline long roundup_over8(long n)
{
    return (n != 0 ? (n < 0? -1 : 1) : 0);
}
static inline long pack8(long n);
static inline long pack8(long n) 
{
    return n < 0 ? n / 8 + roundup_over8(n % 8) : (n + 1) / 8 + roundup_over8((n + 1) % 8) - 1;
}
void packing_test(char *filename, pack_func_t pack_func, unpack_func_t unpack_func, const char *func_name, bool display);




#define NROW(nrl, nrh)  (nrh - nrl + 1)
#define NCOL(ncl, nch)  (nch - ncl + 1)
#endif // __UTIL_H__