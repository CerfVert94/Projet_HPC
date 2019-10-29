#include "mouvement.h"

void routine_FrameDifference(p_image t, p_image t1) {

	int i, j;
	uint8 thresh = THRESHOLD, diff;
	for (i = 0; i < t->height; i++) {
		for (j = 0; j < t->width; j++) {
			diff = abs((t1->img)[i][j].pix - (t->img)[i][j].pix);
			diff > thresh? t->img[i][j].pix = 1 : 0;
		}
	}

}

void SigmaDelta_step1(p_image t, p_image t1) {}

int main(int argc, char const *argv[]){return 0;}

