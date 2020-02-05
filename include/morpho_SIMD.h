/* ------------------------ */
/* --- morpho_SSE2.h --- */
/* ------------------------ */

#ifndef __MORPHO_SSE2_H__
#define __MORPHO_SSE2_H__

#ifndef __MORPHO_SSE2_H__
#define THRESHOLD 20
#define N 2
#endif // __MORPHO_SSE2_H__



#pragma message("  include  morpho_SIMD.h")

void ui8matrix_erosion_SIMD_naive(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);
void ui8matrix_dilation_SIMD_naive(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);

void ui8matrix_erosion_SIMD_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);
void ui8matrix_dilation_SIMD_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);

void ui8matrix_dilation_dilation_SIMD_FO(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);
void ui8matrix_erosion_erosion_SIMD_FO(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);


void ui8matrix_dilation_dilation_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);
void ui8matrix_erosion_erosion_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y);

void test_functions_morpho_SIMD(); 

#endif /* __MORPHO_SSE2_H__ */