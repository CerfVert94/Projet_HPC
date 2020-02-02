/* ---------------- */
/* --- mynrutil.h --- */
/* ---------------- */


#ifndef __MYNRUTIL_H__
#define __MYNRUTIL_H__

#pragma message("  include  mynrutil.h")

#define NR_END 0
#define FREE_ARG char*


/* ------------ */
/* -- memset -- */
/* ------------ */
void memset_ui8matrix(uint8 **X, uint8 value, long nrl, long nrh, long ncl, long nch);

/* ------------ */
/* -- memcmp -- */
/* ------------ */
int memcmp_ui8matrix(uint8 **X, uint8 **Y, long nrl, long nrh, long ncl, long nch);

/* ------------ */
/* -- memcpy -- */
/* ------------ */
void memcpy_ui8matrix(uint8 **X, long nrl, long nrh, long ncl, long nch, uint8 **Y);
uint8** filled_ui8matrix(long nrl, long nrh, long ncl, long nch, uint8 value);

#endif // __NRUTL_H__
