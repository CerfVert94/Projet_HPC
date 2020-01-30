/* ------------------------ */
/* --- morpho_SSE2.h --- */
/* ------------------------ */

#ifndef __MORPHO_SSE2_H__
#define __MORPHO_SSE2_H__

#ifndef __MORPHO_SSE2_H__
#define THRESHOLD 20
#define N 2
#endif // __MORPHO_SSE2_H__



#pragma message("  include  morpho_SSE2.h")

void ui8matrix_erosion_SSE(vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y);
void ui8matrix_dilation_SSE(vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y);

#endif /* __MORPHO_SSE2_H__ */