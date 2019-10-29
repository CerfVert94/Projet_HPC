#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__

#include "img.h"

#define THRESHOLD 5
#define N 2

void routine_FrameDifference(p_image t, p_image t1);

void SigmaDelta_step1(p_image t, p_image t1);

void SigmaDelta_step2(p_image t, p_image t1);

#endif // __MOUVEMENT_H__