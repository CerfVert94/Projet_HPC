/* ------------------------ */
/* --- mouvement_SSE2.c --- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "img_SSE2.h"
#include "mouvement_SSE2.h"
#include "util.h"

/*-----------------------------*/
void SigmaDelta_step0_SSE(p_vimage t0) {
/*-----------------------------*/
	//copy_ui8matrix_ui8matrix(t0->I, t0->nrl, t0->nrh, t0->ncl, t0->nch, t0->M);
	int i, j;
	for (i = 0; i < t0->nrh; i++) {
		for (j = 0; j < t0-> nch; j++)
			t0->V[i][j] = init_vuint8(1);
	}
}

/*-------------------------------------------------*/
void SigmaDelta_step1_SSE(p_vimage vt, p_vimage vt_1) {
/*-------------------------------------------------*/

	vuint8 tI, t_1M, tM, TMP;
	vuint8 C1, C2;

	vuint8 CMP = init_vuint8(127);
	vuint8 ONE = init_vuint8(1);
	vuint8 ZERO = init_vuint8(0);
	vsint8 M_ONE = init_vsint8(-1);

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			tI   = _mm_load_si128((vuint8*) &vt->I[i][j]);
			t_1M = _mm_load_si128((vuint8*) &vt_1->M[i][j]);

			vec_cmplt(t_1M, tI, C1, CMP);
			// tI = _mm_sub_epi8(tI, CMP);
			// t_1M = _mm_sub_epi8(t_1M, CMP);
			// C1  = _mm_cmplt_epi8(t_1M, tI);
			// tI = _mm_add_epi8(tI, CMP);
			// t_1M = _mm_add_epi8(t_1M, CMP);

			C2  = _mm_cmpeq_epi8(tI, t_1M);

		    TMP = _mm_or_si128(_mm_and_si128(C1, ONE), _mm_andnot_si128(C1, M_ONE));
		    TMP = _mm_or_si128(_mm_and_si128(C2, ZERO), _mm_andnot_si128(C2, TMP));
		    tM  = _mm_add_epi8(t_1M, TMP);

		    _mm_store_si128((vuint8*) &vt->M[i][j], tM);

		}
	}
	
}
	

/*----------------------------------*/
void SigmaDelta_step2_SSE(p_vimage vt) {
/*----------------------------------*/

	vuint8 tI, tM, tO;

	vuint8 A, B, C;

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			tI = _mm_load_si128((vuint8*) &vt->I[i][j]);
			tM = _mm_load_si128((vuint8*) &vt->M[i][j]);
			
			tO = vec_subabs(tI, tM);

			_mm_store_si128((vuint8*) &vt->O[i][j], tO);
		}
	}
	
}

/*----------------------------------*/
void SigmaDelta_step3_SSE(p_vimage vt, p_vimage vt_1) {
/*----------------------------------*/

	vuint8 t_1V, tO, tV, TMP;
	vuint8 C1, C2;

	vuint8 CMP = init_vuint8(127);
	vuint8 ONE = init_vuint8(1);
	vuint8 ZERO = init_vuint8(0);
	vsint8 M_ONE = init_vsint8(-1);

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			t_1V = _mm_load_si128((vuint8*) &vt_1->V[i][j]);
			tO 	 = _mm_load_si128((vuint8*) &vt->O[i][j]);
			for (int k = 0; k < N; k++) _mm_add_epi8(tO, tO);

			vec_cmplt(t_1V, tO, C1, CMP);
			// tO = _mm_sub_epi8(tO, CMP);
			// t_1M = _mm_sub_epi8(t_1M, CMP);
			// C1  = _mm_cmplt_epi8(t_1M, tO);
			// tO = _mm_add_epi8(tO, CMP);
			// t_1M = _mm_add_epi8(t_1M, CMP);

			C2  = _mm_cmpeq_epi8(tO, t_1V);

		    TMP = _mm_or_si128(_mm_and_si128(C1, ONE), _mm_andnot_si128(C1, M_ONE));
		    TMP = _mm_or_si128(_mm_and_si128(C2, ZERO), _mm_andnot_si128(C2, TMP));
		    tV  = _mm_add_epi8(t_1V, TMP);

		    _mm_store_si128((vuint8*) &vt->V[i][j], tV);

		}
	}
}

/*----------------------------------*/
void SigmaDelta_step4_SSE(p_vimage vt) {
/*----------------------------------*/

	vuint8 tO, tV, tE;

	vuint8 ZERO = init_vuint8(0);
	vuint8 ONE = init_vuint8(1);

	vuint8 C, CMP = init_vuint8(127);

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			tO = _mm_load_si128((vuint8*) &vt->O[i][j]);
			tV = _mm_load_si128((vuint8*) &vt->V[i][j]);
			
			vec_cmplt(tO, tV, C, CMP);
			// tO = _mm_sub_epi8(tO, CMP);
			// tV = _mm_sub_epi8(tV, CMP);
			// C  = _mm_cmplt_epi8(tO, tV);
			// tO = _mm_add_epi8(tO, CMP);
			// tV = _mm_add_epi8(tV, CMP);

    		tE = _mm_or_si128(_mm_and_si128(C, ZERO), _mm_andnot_si128(C, ONE));
            _mm_store_si128(&vt->E[i][j], tE);
		}
	}
	
}


