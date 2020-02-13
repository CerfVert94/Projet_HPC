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

// void copy_vui8matrix_vui8matrix(vuint8** X, long nrl, long nrh, long vmin, long vmax, vuint8** Y);


/* Init image*/
void SigmaDelta_step0_SIMD(vuint8** I, vuint8** M, vuint8** V, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max); 
void SigmaDelta_step0_SIMD_and_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_step0_SIMD_or_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_step0_SIMD_load_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_step0_SIMD_store_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_step0_SIMD_memset_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_step0_SIMD_InLU_O3_OMP(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_step0_SIMD_InLU_O4_OMP(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);

/* STEP 1 : M Estimation */
void SigmaDelta_step1_SIMD(vuint8** I, vuint8** M_1, vuint8** M, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max); 
/* STEP 2 : O Computation */
void SigmaDelta_step2_SIMD(vuint8** O, vuint8** M, vuint8** I, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max); 
void SigmaDelta_step2_SIMD_SSSE3(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);

void SigmaDelta_step2_InLU_O3_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_step2_ExLU_O3_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max);

/* STEP 3 : V Update & Clamping */
void SigmaDelta_step3_SIMD     (vuint8** V, vuint8** V_1, vuint8** O, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max); 
void SigmaDelta_step3_SIMD_ver2(vuint8** V, vuint8** V_1, vuint8** O, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max); 
void SigmaDelta_step3_InLU_O3_SIMD(vuint8** V, vuint8** V_1, vuint8** O, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max); 
/* STEP 4 : E Estimation */
void SigmaDelta_step4_SIMD(vuint8** O, vuint8** V, vuint8** E, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max); 

void SigmaDelta_SIMD(p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_SIMD_FL(p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_SIMD_FL_OMP(p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
void test_step0_SIMD();

#endif /* __MOUVEMENT_SIMD_H__ */
