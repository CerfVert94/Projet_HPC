/* ---------------- */
/* --- mynrutil.h --- */
/* ---------------- */


#ifndef __MYVNRUTIL_H__
#define __MYVNRUTIL_H__

#pragma message("  include  mynrutil.h")

#define NR_END 0
#define FREE_ARG char*

uint8** vui8matrix_to_ui8matrix(vuint8** vX, long i0, long i1, int j0, int j1, long *nrl, long *nrh, long *ncl, long *nch);

#endif // __NRUTL_H__
