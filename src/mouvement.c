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
	long i, j;
	uint8 thresh = THRESHOLD, diff;
	for (i = t->nrl; i < t->nrh; i++) {
		for (j = t->ncl; j < t->nch; j++) {
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
void SigmaDelta_step0(uint8** M  , uint8** I, uint8** V, long nrl, long nrh, long ncl, long nch) {
/*-----------------------------*/
	copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, M);
	long i, j;
	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++)
			V[i][j] = 1;
	}
}

/*-----------------------------------------*/
void SigmaDelta_step1(uint8** M_1, uint8** I, uint8** M, long nrl, long nrh, long ncl, long nch) {
/*-----------------------------------------*/
	long i, j;

	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++) {
			if (M_1[i][j] < I[i][j])
				M[i][j] = M_1[i][j] + 1;
			else if (M_1[i][j] > I[i][j])
				M[i][j] = M_1[i][j] - 1;
			else
				M[i][j] = M_1[i][j];
		}
	}

}

/*-----------------------------------------*/
void SigmaDelta_step2(uint8** M  , uint8** I, uint8** O, long nrl, long nrh, long ncl, long nch) {
/*-----------------------------------------*/
	long i, j;

	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++)
			O[i][j] = abs(M[i][j] - I[i][j]);
	}

}

/*-----------------------------------------*/
void SigmaDelta_step3(uint8** V_1, uint8** O, uint8** V, long nrl, long nrh, long ncl, long nch) {
/*-----------------------------------------*/
	long i, j;

	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++) {
			if (V_1[i][j] < (N * O[i][j]))
				V[i][j] = V_1[i][j] + 1;
			else if (V_1[i][j] > (N * O[i][j]))
				V[i][j] = V_1[i][j] - 1;
			else
				V[i][j] = V_1[i][j];
			V[i][j] = max(min(V[i][j], Vmax), Vmin);
		}
	}

}

/*-----------------------------------------*/
void SigmaDelta_step4(uint8** O, uint8** V, uint8** E, long nrl, long nrh, long ncl, long nch) {
/*-----------------------------------------*/
	long i, j;

	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++) {
			if (O[i][j] < V[i][j] )
				E[i][j] = 0;
			else
				E[i][j] = 1;
		}
	}

}

/*-----------------------------------------*/
void SigmaDelta(p_image t, p_image t_1) {
/*-----------------------------------------*/
	SigmaDelta_step1(t->I, t_1->M, t->M, t->nrl, t->nrh, t->ncl, t->nch);
	SigmaDelta_step2(t->O, t->M, t->I, t->nrl, t->nrh, t->ncl, t->nch);
	SigmaDelta_step3(t->V, t_1->V, t->O, t->nrl, t->nrh, t->ncl, t->nch);
	SigmaDelta_step4(t->O, t->V, t->E, t->nrl, t->nrh, t->ncl, t->nch);

}

/*-------*/
void test_mouvement() {
/*-------*/
	
	p_image t_1 = create_image("../car3/car_3000.pgm");
	p_image t = create_image("../car3/car_3001.pgm");

	printf("Nrh: %ld\n", t->nrh);
	printf("Nch: %ld\n", t->nch);

	SigmaDelta_step0(t_1->I, t_1->M, t_1->V, t->nrl, t->nrh, t->ncl, t->nch);
	SigmaDelta(t, t_1);
	for (long i = t->nrl; i < 50; i++) { 
		for (long j = t->ncl; j < 50; j++) {
			printf("%d ", t->E[i][j]);
		}
		printf("\n");
	}

	free_image(t);
	free_image(t_1);
}