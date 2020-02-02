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

/*-----------------------------*/
void SigmaDelta_step0_mem(uint8** M  , uint8** I, uint8** V, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------*/
	memcpy_ui8matrix(I, nrl, nrh, ncl, nch, M);
	memset_ui8matrix(V, v_min, nrl, nrh, ncl, nch);
}

#define ESTIMATE_M(M_1, I, M, j)\
M[(j)] = M_1[(j)] + (M_1[(j)] < I[(j)]) - (M_1[(j)] > I[(j)]); \
/*-----------------------------------------*/
void SigmaDelta_step1_InLU_O3_NoIf(uint8** M_1, uint8** I, uint8** M, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j, r;
	const long order = 3;
	uint8 *M_1_row0, *M_1_row1, *M_1_row2;
	uint8 *M_row0  , *M_row1  , *M_row2;
	uint8 *I_row0  , *I_row1  , *I_row2;
	r = (nch + 1) % order;

	for (i = nrl; i <= nrh; i++) {
		M_1_row0 = M_1[i + 0];I_row0 = I[i + 0]; M_row0 = M[i + 0];
		for (j = ncl; j <= nch - r; j += order) {
			ESTIMATE_M(M_1_row0, I_row0, M_row0, j + 0);
			ESTIMATE_M(M_1_row0, I_row0, M_row0, j + 1);
			ESTIMATE_M(M_1_row0, I_row0, M_row0, j + 2);
		}
	}
	switch (r) {
		case 2:
		for (i = nrl; i <= nrh; i++) {
			M_1_row0 = M_1[i + 0];I_row0 = I[i + 0]; M_row0 = M[i + 0];
			ESTIMATE_M(M_1_row0, I_row0, M_row0, nch - 1);
			ESTIMATE_M(M_1_row0, I_row0, M_row0, nch + 0);
		}
		break;
		case 1:
		for (i = nrl; i <= nrh; i++) {
			M_1_row0 = M_1[i + 0];I_row0 = I[i + 0]; M_row0 = M[i + 0];
			ESTIMATE_M(M_1_row0, I_row0, M_row0, nch + 0);
		}
		break;
		case 0:
		break;
	}
}

void SigmaDelta_step1_ExLU_O3_NoIf(uint8** M_1, uint8** I, uint8** M, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j, r;
	const long order = 3;
	uint8 *M_1_row0, *M_1_row1, *M_1_row2;
	uint8 *M_row0  , *M_row1  , *M_row2;
	uint8 *I_row0  , *I_row1  , *I_row2;
	r = (nrh + 1) % order;

	for (i = nrl; i <= nrh - r; i += order) {
		M_1_row0 = M_1[i + 0]; I_row0 = I[i + 0]; M_row0 = M[i + 0];
		M_1_row1 = M_1[i + 1]; I_row1 = I[i + 1]; M_row1 = M[i + 1];
		M_1_row2 = M_1[i + 2]; I_row2 = I[i + 2]; M_row2 = M[i + 2];
		for (j = ncl; j <= nch; j ++) {
			ESTIMATE_M(M_1_row0, I_row0, M_row0, j);
			ESTIMATE_M(M_1_row1, I_row1, M_row1, j);
			ESTIMATE_M(M_1_row2, I_row2, M_row2, j);
		}
	}
	switch (r) {
		case 2:
			M_1_row0 = M_1[nrh - 1]; I_row0 = I[nrh - 1]; M_row0 = M[nrh - 1];
			M_1_row1 = M_1[nrh - 0]; I_row1 = I[nrh - 0]; M_row1 = M[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				ESTIMATE_M(M_1_row0, I_row0, M_row0, j);
				ESTIMATE_M(M_1_row1, I_row1, M_row1, j);
			}
			break;
		case 1:
			M_1_row1 = M_1[nrh - 0]; I_row1 = I[nrh - 0]; M_row1 = M[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				ESTIMATE_M(M_1_row1, I_row1, M_row1, j);
			}
			break;
		case 0:
		break;
	}
}

/*-----------------------------------------*/
void SigmaDelta_step2_InLU_O3_bitop(uint8** M  , uint8** I, uint8** O, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j, r;
	const long order = 3;
	uint8 *M_row0  , *M_row1  , *M_row2;
	uint8 *O_row0  , *O_row1  , *O_row2;
	uint8 *I_row0  , *I_row1  , *I_row2;
	int16 x0, x1, x2, y0, y1, y2;
	r = (nch + 1) % order;

	for (i = nrl; i <= nrh; i++) {
		M_row0 = M[i + 0];I_row0 = I[i + 0]; O_row0 = O[i + 0];
		for (j = ncl; j <= nch - r; j += order) {
			x0 = M_row0[j + 0] - I_row0[j + 0];
			x1 = M_row0[j + 1] - I_row0[j + 1];
			x2 = M_row0[j + 2] - I_row0[j + 2];
			y0 = x0 >> 8;
			y1 = x1 >> 8;
			y2 = x2 >> 8;
			O_row0[j + 0] = (y0 + x0) ^ y0;
			O_row0[j + 1] = (y1 + x1) ^ y1;
			O_row0[j + 2] = (y2 + x2) ^ y2;
		}
	}
	switch (r) {
		case 2:
		for (i = nrl; i <= nrh; i++) {
			M_row0 = M[i + 0];I_row0 = I[i + 0]; O_row0 = O[i + 0];
			x0 = M_row0[nch - 1] - I_row0[nch - 1];
			x1 = M_row0[nch - 0] - I_row0[nch - 0];
			y0 = x0 >> 8;
			y1 = x1 >> 8;
			O_row0[nch - 1] = (y0 + x0) ^ y0;
			O_row0[nch - 0] = (y1 + x1) ^ y1;
		}
		break;
		case 1:
		
		for (i = nrl; i <= nrh; i++) {
			x1 = M_row0[nch - 0] - I_row0[nch - 0];
			y1 = x1 >> 8;
			O_row0[nch - 0] = (y1 + x1) ^ y1;
		}
		break;
		case 0:
		break;
	}
}
void SigmaDelta_step2_ExLU_O3_bitop(uint8** M  , uint8** I, uint8** O, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {

	long i, j, r;
	const long order = 3;
	uint8 *M_row0  , *M_row1  , *M_row2;
	uint8 *O_row0  , *O_row1  , *O_row2;
	uint8 *I_row0  , *I_row1  , *I_row2;
	int16 x0, x1, x2, y0, y1, y2;
	
	r = (nrh + 1) % order;

	for (i = nrl; i <= nrh - r; i += order) {
		M_row0 = M[i + 0];I_row0 = I[i + 0]; O_row0 = O[i + 0];
		M_row1 = M[i + 1];I_row1 = I[i + 1]; O_row1 = O[i + 1];
		M_row2 = M[i + 2];I_row2 = I[i + 2]; O_row2 = O[i + 2];
		for (j = ncl; j <= nch; j ++) {
			x0 = M_row0[j] - I_row0[j];
			x1 = M_row1[j] - I_row1[j];
			x2 = M_row2[j] - I_row2[j];
			y0 = x0 >> 8;
			y1 = x1 >> 8;
			y2 = x2 >> 8;
			O_row0[j] = (y0 + x0) ^ y0;
			O_row1[j] = (y1 + x1) ^ y1;
			O_row2[j] = (y2 + x2) ^ y2;
		}
	}
	switch (r) {
		case 2:
			O_row0 = O[nrh - 1]; I_row0 = I[nrh - 1]; M_row0 = M[nrh - 1];
			O_row1 = O[nrh - 0]; I_row1 = I[nrh - 0]; M_row1 = M[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				x0 = M_row0[j] - I_row0[j];
				x1 = M_row1[j] - I_row1[j];
				y0 = x0 >> 8;
				y1 = x1 >> 8;
				O_row0[j] = (y0 + x0) ^ y0;
				O_row1[j] = (y1 + x1) ^ y1;
			}
			break;
		case 1:
			O_row1 = O[nrh - 0]; I_row1 = I[nrh - 0]; M_row1 = M[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				x1 = M_row1[j] - I_row1[j];
				y1 = x1 >> 8;
				O_row1[j] = (y1 + x1) ^ y1;
			}
			break;
		case 0:
		break;
	}
}

#define UPDATE_CLAMP_V(V_1, N, O, V, j, v_t,v_min, v_max)\
v_t = V_1[(j)];\
v_t += V_1[(j)] < N * O[(j)];\
v_t -= V_1[(j)] > N * O[(j)];\
v_t = (v_t > v_max ? v_max : v_t);\
v_t = (v_t < v_min ? v_min : v_t);\
V[(j)] = v_t;					//max(min(v_t, v_max), v_min);


/*-----------------------------------------*/
void SigmaDelta_step3_InLU_O3_NoIf(uint8** V_1, uint8** O, uint8** V, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j, r;
	const long order = 3;
	uint8 *V_1_row0, *V_1_row1, *V_1_row2;
	uint8 *V_row0  , *V_row1  , *V_row2;
	uint8 *O_row0  , *O_row1  , *O_row2;
	int16 v_t;
	r = (nch + 1) % order;

	for (i = nrl; i <= nrh; i++) {
		V_1_row0 = V_1[i + 0];O_row0 = O[i + 0]; V_row0 = V[i + 0];
		for (j = ncl; j <= nch - r; j += order) {
			UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, j + 0, v_t, v_min, v_max);
			UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, j + 1, v_t, v_min, v_max);
			UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, j + 2, v_t, v_min, v_max);
		}
	}
	switch (r) {
		case 2:
		for (i = nrl; i <= nrh; i++) {
			V_1_row0 = V_1[i + 0];O_row0 = O[i + 0]; V_row0 = V[i + 0];
			UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, nch - 1, v_t, v_min, v_max);
			UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, nch + 0, v_t, v_min, v_max);
		}
		break;
		case 1:
		for (i = nrl; i <= nrh; i++) {
			V_1_row0 = V_1[i + 0];O_row0 = O[i + 0]; V_row0 = V[i + 0];
			UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, nch + 0, v_t, v_min, v_max);
		}
		break;
		case 0:
		break;
	}
}

/*-----------------------------------------*/
void SigmaDelta_step3_ExLU_O3_NoIf(uint8** V_1, uint8** O, uint8** V, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j, r;
	const long order = 3;
	uint8 *V_1_row0, *V_1_row1, *V_1_row2;
	uint8 *V_row0  , *V_row1  , *V_row2;
	uint8 *O_row0  , *O_row1  , *O_row2;
	int16 v_t;
	r = (nrh + 1) % order;

	for (i = nrl; i <= nrh - r; i += order) {
		V_1_row0 = V_1[i + 0]; O_row0 = O[i + 0]; V_row0 = V[i + 0];
		V_1_row1 = V_1[i + 1]; O_row1 = O[i + 1]; V_row1 = V[i + 1];
		V_1_row2 = V_1[i + 2]; O_row2 = O[i + 2]; V_row2 = V[i + 2];
		for (j = ncl; j <= nch; j ++) {
			UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, j + 0, v_t, v_min, v_max);
			UPDATE_CLAMP_V(V_1_row1, n_coeff, O_row1, V_row1, j + 0, v_t, v_min, v_max);
			UPDATE_CLAMP_V(V_1_row2, n_coeff, O_row2, V_row2, j + 0, v_t, v_min, v_max);
		}
	}
	switch (r) {
		case 2:
			V_1_row0 = V_1[nrh - 1]; O_row0 = O[nrh - 1]; V_row0 = V[nrh - 1];
			V_1_row1 = V_1[nrh - 0]; O_row1 = O[nrh - 0]; V_row1 = V[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				UPDATE_CLAMP_V(V_1_row0, n_coeff, O_row0, V_row0, j + 0, v_t, v_min, v_max);
				UPDATE_CLAMP_V(V_1_row1, n_coeff, O_row1, V_row1, j + 0, v_t, v_min, v_max);
			}
			break;
		case 1:
			V_1_row1 = V_1[nrh - 0]; O_row1 = O[nrh - 0]; V_row1 = V[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				UPDATE_CLAMP_V(V_1_row1, n_coeff, O_row1, V_row1, j + 0, v_t, v_min, v_max);
			}
			break;
		case 0:
		break;
	}
}

/*-----------------------------------------*/
void SigmaDelta_step4_InLU_O3_NoIf(uint8** O, uint8** V, uint8** E, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j, r;
	const long order = 3;
	uint8 *E_row0  , *E_row1  , *E_row2;
	uint8 *V_row0  , *V_row1  , *V_row2;
	uint8 *O_row0  , *O_row1  , *O_row2;
	int16 v_t;
	r = (nch + 1) % order;

	for (i = nrl; i <= nrh; i++) {
		E_row0 = E[i + 0]; O_row0 = O[i + 0]; V_row0 = V[i + 0];
		for (j = ncl; j <= nch - r; j += order) {
			E_row0[j + 0] = O_row0[j + 0] >= V_row0[j + 0];
			E_row0[j + 1] = O_row0[j + 1] >= V_row0[j + 1];
			E_row0[j + 2] = O_row0[j + 2] >= V_row0[j + 2];
		}
	}
	switch (r) {
		case 2:
		for (i = nrl; i <= nrh; i++) {
			E_row0 = E[i + 0];O_row0 = O[i + 0]; V_row0 = V[i + 0];
			E_row0[nch + 0] = O_row0[nch + 0] >= V_row0[j + 0];
			E_row0[nch - 1] = O_row0[nch - 1] >= V_row0[j - 1];
		}
		break;
		case 1:
		for (i = nrl; i <= nrh; i++) {
			E_row0 = E[i + 0];O_row0 = O[i + 0]; V_row0 = V[i + 0];
			E_row0[nch + 0] = O_row0[nch + 0] >= V_row0[j + 0];
		}
		break;
		case 0:
		break;
	}
}

void SigmaDelta_step4_ExLU_O3_NoIf(uint8** O, uint8** V, uint8** E, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------*/
	long i, j, r;
	const long order = 3;
	uint8 *E_row0  , *E_row1  , *E_row2;
	uint8 *V_row0  , *V_row1  , *V_row2;
	uint8 *O_row0  , *O_row1  , *O_row2;
	int16 v_t;
	r = (nrh + 1) % order;

	// memset_ui8matrix(E, 1, nrl, nrh, ncl, nch);
	for (i = nrl; i <= nrh - r; i += order) {
		E_row0 = E[i + 0]; O_row0 = O[i + 0]; V_row0 = V[i + 0];
		E_row1 = E[i + 1]; O_row1 = O[i + 1]; V_row1 = V[i + 1];
		E_row2 = E[i + 2]; O_row2 = O[i + 2]; V_row2 = V[i + 2];
		for (j = ncl; j <= nch; j ++) {
				E_row0[j] = O_row0[j] >= V_row0[j];
				E_row1[j] = O_row1[j] >= V_row1[j];
				E_row2[j] = O_row2[j] >= V_row2[j];
		}
	}
	switch (r) {
		case 2:
			E_row0 = E[nrh - 1]; O_row0 = O[nrh - 1]; V_row0 = V[nrh - 1];
			E_row1 = E[nrh - 0]; O_row1 = O[nrh - 0]; V_row1 = V[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				E_row0[j] = O_row0[j] >= V_row0[j];
				E_row1[j] = O_row1[j] >= V_row1[j];
			}
			break;
		case 1:
			E_row1 = E[nrh - 0]; O_row1 = O[nrh - 0]; V_row1 = V[nrh - 0];
			for (j = ncl; j <= nch; j++) {
				E_row1[j] = O_row1[j] >= V_row1[j];
			}
			break;
		case 0:
		break;
	}
}


void SigmaDelta_best(p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	SigmaDelta_step1_InLU_O3_NoIf (t0->M, t1->I, t1->M, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
	SigmaDelta_step2_InLU_O3_bitop(t1->M, t1->I, t1->O, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
	SigmaDelta_step3_InLU_O3_NoIf (t0->V, t1->O, t1->V, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
	SigmaDelta_step4_InLU_O3_NoIf (t1->O, t1->V, t1->E, t1->nrl, t1->nrh, t1->ncl, t1->nch, n_coeff, v_min, v_max);
}
