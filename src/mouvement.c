/* ------------------------ */
/* ----- mouvement.c ------ */
/* ------------------------ */

#include <stdlib.h>
#include <stdio.h>

#include "nrdef.h"
#include "nrutil.h"

#include "util.h"
#include "img.h"

#include "mouvement.h"


/*-----------------------------------------------*/
void routine_FrameDifference(p_image t, p_image t1) {
/*-----------------------------------------------*/
	int i, j;
	uint8 thresh = THRESHOLD, diff;
	for (i = 0; i < t->nrh; i++) {
		for (j = 0; j < t->nch; j++) {
			diff = abs(t1->I[i][j] - t->I[i][j]);
			if (diff > thresh) {
				t->E[i][j] = 1;
			}
			else {
				t->E[i][j] = 0;
			}
		}
	}

}

/*-----------------------------*/
void SigmaDelta_step0(p_image t0) {
/*-----------------------------*/
	copy_ui8matrix_ui8matrix(t0->I, t0->nrl, t0->nrh, t0->ncl, t0->nch, t0->M);
	int i, j;
	for (i = 0; i < t0->nrh; i++) {
		for (j = 0; j < t0-> nch; j++)
			t0->V[i][j] = 1;
	}
}

/*-----------------------------------------*/
void SigmaDelta_step1(p_image t, p_image t_1) {
/*-----------------------------------------*/
	int i, j;

	for (i = 0; i < t->nrh; i++) {
		for (j = 0; j < t->nch; j++) {
			if (t_1->M[i][j] < t->I[i][j])
				t->M[i][j] = t_1->M[i][j] + 1;
			else if (t_1->M[i][j] > t->I[i][j])
				t->M[i][j] = t_1->M[i][j] - 1;
			else
				t->M[i][j] = t_1->M[i][j];
		}
	}

}

/*-----------------------------------------*/
void SigmaDelta_step2(p_image t) {
/*-----------------------------------------*/
	int i, j;

	for (i = 0; i < t->nrh; i++) {
		for (j = 0; j < t->nch; j++)
			t->O[i][j] = abs(t->M[i][j] - t->I[i][j]);
	}

}

/*-----------------------------------------*/
void SigmaDelta_step3(p_image t, p_image t_1) {
/*-----------------------------------------*/
	int i, j;

	for (i = 0; i < t->nrh; i++) {
		for (j = 0; j < t->nch; j++) {
			if (t_1->V[i][j] < (N * t->O[i][j]))
				t->V[i][j] = t_1->V[i][j] + 1;
			else if (t_1->V[i][j] > (N * t->O[i][j]))
				t->V[i][j] = t_1->V[i][j] - 1;
			else
				t->V[i][j] = t_1->V[i][j];
			t->V[i][j] = max(min(t->V[i][j], Vmax), Vmin);
		}
	}

}

/*-----------------------------------------*/
void SigmaDelta_step4(p_image t) {
/*-----------------------------------------*/
	int i, j;

	for (i = 0; i < t->nrh; i++) {
		for (j = 0; j < t->nch; j++) {
			if (t->O[i][j] < t->V[i][j] )
				t->E[i][j] = 0;
			else
				t->E[i][j] = 1;
		}
	}

}


/*-----------------------------------------*/
void SigmaDelta(p_image t, p_image t_1) {
/*-----------------------------------------*/
	SigmaDelta_step1(t, t_1);
	SigmaDelta_step2(t);
	SigmaDelta_step3(t, t_1);
	SigmaDelta_step4(t);

}


/*-------*/
void test() {
/*-------*/
	
	p_image t_1 = create_image("../car3/car_3000.pgm");
	p_image t = create_image("../car3/car_3001.pgm");

	printf("Nrh: %ld\n", t->nrh);
	printf("Nch: %ld\n", t->nch);

	SigmaDelta_step0(t_1);
	SigmaDelta(t, t_1);
	for (int i = t->nrl; i < 50; i++) { 
		for (int j = t->ncl; j < 50; j++) {
			printf("%d ", t->E[i][j]);
		}
		printf("\n");
	}

	free_image(t);
	free_image(t_1);
}