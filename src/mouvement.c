/* ------------------------ */
/* ----- mouvement.c ------ */
/* ------------------------ */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include "nrdef.h"
#include "vnrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
#include "vnrutil.h"
#include "myvnrutil.h"

#include "util.h"
#include "img.h"
#include "img_SIMD.h"


#include "mouvement.h"


/*-----------------------------------------------*/
void routine_FrameDifference(p_image t0, p_image t1) {
/*-----------------------------------------------*/
	long i, j;
	uint8 thresh = THRESHOLD, diff;
	for (i = t0->nrl; i < t0->nrh; i++) {
		for (j = t0->ncl; j < t0->nch; j++) {
			diff = abs(t1->I[i][j] - t0->I[i][j]);
			if (diff > thresh) {
				t0->E[i][j] = 1;
			}
			else {
				t0->E[i][j] = 0;
			}
		}
	}

}

/*-----------------------------*/
void SigmaDelta_step0_naive(uint8** M  , uint8** I, uint8** V, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------*/
	copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, M);
	// display_ui8matrix(M, nrl, nrh, ncl, nch, "%u", "M_NAIVE");
	long i, j;
	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++)
			V[i][j] = v_min;
	}
}

/*-----------------------------------------*/
void SigmaDelta_step1_naive(uint8** M_1, uint8** I, uint8** M, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j;

	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++) {
			if ((M_1[i][j] < I[i][j]) && (M_1[i][j] <= Vmax))
				M[i][j] = M_1[i][j] + 1;
			else if ((M_1[i][j] > I[i][j]) && (M_1[i][j] >= Vmin))
				M[i][j] = M_1[i][j] - 1;
			else
				M[i][j] = M_1[i][j];
		}
	}
}


/*-----------------------------------------*/
void SigmaDelta_step2_naive(uint8** M  , uint8** I, uint8** O, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j;

	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++)
			O[i][j] = abs(M[i][j] - I[i][j]);
	}

}
/*-----------------------------------------*/
void SigmaDelta_step3_naive(uint8** V_1, uint8** O, uint8** V, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j;
	uint8 v_0;
	for (i = nrl; i <= nrh; i++) {
		for (j = ncl; j <= nch; j++) {
			
			if ((V_1[i][j] < (n_coeff * O[i][j])) && V_1[i][j] <= Vmax)
				v_0 = V_1[i][j] + 1;
			else if ((V_1[i][j] > (n_coeff * O[i][j])) && V_1[i][j] >= Vmin)
				v_0 = V_1[i][j] - 1;
			else
				v_0 = V_1[i][j];
			V[i][j] = max(min(v_0, v_max), v_min);
		}
	}

}
/*-----------------------------------------*/
void SigmaDelta_step4_naive(uint8** O, uint8** V, uint8** E, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
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
void SigmaDelta_naive(p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	SigmaDelta_step1_naive(t0->M, t1->I, t1->M, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
	SigmaDelta_step2_naive(t1->M, t1->I, t1->O, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
	SigmaDelta_step3_naive(t0->V, t1->O, t1->V, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
	SigmaDelta_step4_naive(t1->O, t1->V, t1->E, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
}
