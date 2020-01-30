/* ------------------------ */
/* ----- mouvement.h ------ */
/* ------------------------ */

#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__

#define THRESHOLD 20
#define N 2

void routine_FrameDifference(p_image t, p_image t1);

/* Init image*/
void SigmaDelta_step0(uint8** I, uint8** M, uint8** V, long nrl, long nrh, long ncl, long nch);
/* STEP 1 : M Estimation */
void SigmaDelta_step1(uint8** I, uint8** M_1, uint8** M, long nrl, long nrh, long ncl, long nch);
/* STEP 2 : O Computation */
void SigmaDelta_step2(uint8** O, uint8** M, uint8** I, long nrl, long nrh, long ncl, long nch);
/* STEP 3 : V Update & Clamping */
void SigmaDelta_step3(uint8** V, uint8** V_1, uint8** O, long nrl, long nrh, long ncl, long nch);
/* STEP 4 : E Estimation */
void SigmaDelta_step4(uint8** O, uint8** V, uint8** E, long nrl, long nrh, long ncl, long nch);
void SigmaDelta(p_image t, p_image t_1);



void test_mouvement();

#endif // __MOUVEMENT_H__