/* ------------------------ */
/* --- mouvement_SIMD.h --- */
/* ------------------------ */

#ifndef __MOUVEMENT_SIMD_H__
#define __MOUVEMENT_SIMD_H__

#ifndef __MOUVEMENT_H__
#define THRESHOLD 20
#define N 2
#endif // __MOUVEMENT_H__



#pragma message("  include  mouvement_SIMD.h")

void copy_vui8matrix_vui8matrix(vuint8** X, long nrl, long nrh, long vmin, long vmax, vuint8** Y);


/* Init image*/
void SigmaDelta_step0_SIMD(vuint8** I, vuint8** M, vuint8** V, long nrl, long nrh, int v0, int v1); 
/* STEP 1 : M Estimation */
void SigmaDelta_step1_SIMD(vuint8** I, vuint8** M_1, vuint8** M, long nrl, long nrh, int v0, int v1); 
/* STEP 2 : O Computation */
void SigmaDelta_step2_SIMD(vuint8** O, vuint8** M, vuint8** I, long nrl, long nrh, int v0, int v1); 
/* STEP 3 : V Update & Clamping */
void SigmaDelta_step3_SIMD(vuint8** V, vuint8** V_1, vuint8** O, long nrl, long nrh, int v0, int v1); 
/* STEP 4 : E Estimation */
void SigmaDelta_step4_SIMD(vuint8** O, vuint8** V, vuint8** E, long nrl, long nrh, int v0, int v1); 

void test_step0_SIMD();

#endif /* __MOUVEMENT_SIMD_H__ */
