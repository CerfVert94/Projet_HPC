#include <stdlib.h>
#include <stdio.h>
#include "nrdef.h"
#include "img.h"
#include "mouvement.h"

void routine_FrameDifference(p_image t, p_image t1) {

	int i, j;
	uint8 thresh = THRESHOLD, diff;
	for (i = 0; i < t->nrl; i++) {
		for (j = 0; j < t->ncl; j++) {
			diff = abs((t1->img)[i][j] - (t->img)[i][j]);
			if (diff > thresh)
				t->move[i][j] = 1;
			else
				t->move[i][j] = 0;
		}
	}

}

void SigmaDelta_step1(p_image t, p_image t1) {}

void test() {
	p_image t = create_image("../car3/car_3000.pgm");
	p_image t1 = create_image("../car3/car_3001.pgm");

	printf("Nrh: %ld\n", t->nrh);
	printf("Nch: %ld\n", t->nch);
	routine_FrameDifference(t, t1);
	// for (int i = 150; i < 200; i++) { 
	// 	for (int j = 150; j < 200; j++) {
	// 		printf("%d ", t->move[i][j]);
	// 	}
	// 	printf("\n");
	// }

	free_image(t);
	free_image(t1);
}