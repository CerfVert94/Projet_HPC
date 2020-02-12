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
void ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR_OMP(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z);
void ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR_OMP(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z);
void ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z);
void ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z);
void ui8matrix_sequence_SIMD_Pipeline2_FO_InLU_O3_ValAddrRR(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z);

void ui8matrix_erosion_SIMD_divide_row_and_conquer(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);
void ui8matrix_erosion_SIMD_divide_col_and_conquer(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);
void ui8matrix_erosion_SIMD_pipeline2_LU3x3_InLU_O3_RR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);
void ui8matrix_erosion_SIMD_InLU_O3_AddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) ;
void ui8matrix_erosion_SIMD_col_pipeline(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);

void ui8matrix_erosion_SIMD_naive(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);
void ui8matrix_dilation_SIMD_naive(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);

void ui8matrix_erosion_SIMD_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);
void ui8matrix_dilation_SIMD_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);

void ui8matrix_dilation_dilation_SIMD_FO(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);
void ui8matrix_erosion_erosion_SIMD_FO(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);

void ui8matrix_dilation_SIMD_InLU_O3_AddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);
void ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);
void ui8matrix_dilation_SIMD_InLU_O3_ValAddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);
void ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);

void ui8matrix_dilation_dilation_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);
void ui8matrix_erosion_erosion_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer, vuint8 **Y);


void ui8matrix_erosion_SIMD_InLU_O3_AddrRR_OMP (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y);
void test_functions_morpho_SIMD(); 

#endif /* __MORPHO_SSE2_H__ */