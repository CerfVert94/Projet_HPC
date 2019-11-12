#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__


#define THRESHOLD 20
#define N 3

void routine_FrameDifference(p_image t, p_image t1);


void SigmaDelta_step0(p_image t0); /* Init image at time 0 */
void SigmaDelta_step1(p_image t, p_image t_1); /* STEP 1 : M Estimation */
void SigmaDelta_step2(p_image t); /* STEP 2 : O Computation */
void SigmaDelta_step3(p_image t, p_image t_1); /* STEP 3 : V Update & Clamping */
void SigmaDelta_step4(p_image t); /* STEP 4 : E Estimation */
void SigmaDelta();

void test();

#endif // __MOUVEMENT_H__