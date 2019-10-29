#include "mouvement.h"

void routine_FrameDifference(p_image t, p_image t1) {

	int i, j;
	uint8 tresh = THRESHOLD, diff;
	for (i = 0; i < t->height; i++) {
		for (j = 0; j < t->weight; j++) {
			diff = abs((t1->img)[i][j] - (t->img)[i][j]);
			diff > thresh? t->img[i][j] = 1 : ;
		}
	}

}

void SigmaDelta_step1(p_image t, p_image t1) {}

int main(int argc, char const *argv[]){return 0;}

