/* ---------------- */
/* --- mynrutil.h --- */
/* ---------------- */


#ifndef __MYVNRUTIL_H__
#define __MYVNRUTIL_H__

#pragma message("  include  mynrutil.h")

#define NR_END 0
#define FREE_ARG char*

uint8  **vui8matrix_to_ui8matrix(vuint8 **vX, long i0, long i1, int j0, int j1, long *nrl, long *nrh, long *ncl, long *nch);
vuint8 **ui8matrix_to_vui8matrix(uint8  **X , long  nrl, long  nrh, long ncl, long nch, int *i0, int *i1, int *j0, int *j1);

vuint8 **LoadPGM_vui8matrix(char *filename, long *nrl, long *nrh, int *v0, int *v1);
#endif // __NRUTL_H__
