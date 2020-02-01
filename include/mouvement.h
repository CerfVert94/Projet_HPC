/* ------------------------ */
/* ----- mouvement.h ------ */
/* ------------------------ */

#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__

#define THRESHOLD 20
#define N 2

struct sd_set{
    char func_name[128];
    void (*sd_func)(uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
    uint8 n_coeff; 
    uint8 v_min; 
    uint8 v_max; 
    // enum {SD_STEP0, SD_STEP1, SD_STEP2, SD_STEP3, SD_STEP4}sd_type;
};
struct complete_sd_set{
    char func_name[128];
    void (*sd_step0)(uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
    void (*sd_func)(p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
    uint8 n_coeff;
    uint8 v_min; 
    uint8 v_max; 
};



void routine_FrameDifference(p_image t, p_image t1);

/* Init image*/
void SigmaDelta_step0_naive(uint8** M  , uint8** I, uint8** V, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
/* STEP 1 : M Estimation */
void SigmaDelta_step1_naive(uint8** M_1, uint8** I, uint8** M, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
/* STEP 2 : O Computation */
void SigmaDelta_step2_naive(uint8** M  , uint8** I, uint8** O, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
/* STEP 3 : V Update & Clamping */
void SigmaDelta_step3_naive(uint8** V_1, uint8** O, uint8** V, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
/* STEP 4 : E Estimation */
void SigmaDelta_step4_naive(uint8** O, uint8** V, uint8** E, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
void SigmaDelta_naive(p_image t, p_image t_1, uint8 n_coeff, uint8 v_min, uint8 v_max);



void test_mouvement();

#endif // __MOUVEMENT_H__