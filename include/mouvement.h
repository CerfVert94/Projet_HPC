#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__


#define THRESHOLD 20
#define N 2

void routine_FrameDifference(p_image t, p_image t1);

void init_SigmaDelta(p_image t);

void SigmaDelta_step0(p_image t, p_image t1);

void SigmaDelta_step1(p_image t, p_image t1);

void test();

#endif // __MOUVEMENT_H__