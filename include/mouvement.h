#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__


#define THRESHOLD 20
#define N 3

void routine_FrameDifference(p_image t, p_image t1);

/* Init image at time 0 */
void SigmaDelta_step0(p_image t0);

/* SigmaDelta algorithm between image t_1 and next time image t */
void SigmaDelta_step1(p_image t, p_image t_1);

void test();

#endif // __MOUVEMENT_H__