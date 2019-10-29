#include "nrdef.h"
#include "img.h"
#include "mouvement.h"
#include <stdlib.h>

void routine_FrameDifference(p_image t, p_image t1) {

	int i, j;
	uint8 thresh = THRESHOLD, diff;
	for (i = 0; i < t->height; i++) {
		for (j = 0; j < t->width; j++) {
			diff = abs(&(t1->img)[i][j] - &(t->img)[i][j]);
			if (diff > thresh) 
				 t->img[i][j].move = 1;
		}
	}

}

void SigmaDelta_step1(p_image t, p_image t1) {}

//int main(int argc, char const *argv[]){return 0;}

