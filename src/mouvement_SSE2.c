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

/*-------------------------------------------------*/
void SigmaDelta_step1_SSE(p_vimage vt, p_vimage vt_1) {
/*-------------------------------------------------*/

	vuint8 tI, t_1I;
	vuint8 c, S;

	vuint8 SUB = init_vuint8(127);

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			tI   = _mm_load_si128((vuint8*) &vt->I[i][j]);
			t_1I = _mm_load_si128((vuint8*) &vt_1->I[i][j]);

			vec_cmplsb(tI, t_1I, c, SUB);
			


		}
	}
	
}
	

/*----------------------------------*/
void SigmaDelta_step2_SSE(p_vimage vt) {
/*----------------------------------*/

	vuint8 tI, tM;
	vuint8 vsub;

	vuint8 A, B, C;

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			tI = _mm_load_si128((vuint8*) &vt->I[i][j]);
			tM = _mm_load_si128((vuint8*) &vt->M[i][j]);
			
			vsub = vec_subabs(tI, tM);

			_mm_store_si128((vuint8*) &vt->O[i][j], vsub);
		}
	}
	
}

/*----------------------------------*/
void SigmaDelta_step4_SSE(p_vimage vt) {
/*----------------------------------*/

	vuint8 tO, tV, tE;

	vuint8 ZERO = init_vuint8(0);
	vuint8 ONE = init_vuint8(1);

	vuint8 C, MSB = init_vuint8(127);

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			tO = _mm_load_si128((vuint8*) &vt->O[i][j]);
			tV = _mm_load_si128((vuint8*) &vt->V[i][j]);
			
			tO = _mm_sub_epi8(tO, MSB);
			tV = _mm_sub_epi8(tV, MSB);
			C  = _mm_cmplt_epi8(tO, tV);
			tO = _mm_add_epi8(tO, MSB);
			tV = _mm_add_epi8(tV, MSB);

    		tE = _mm_or_si128(_mm_and_si128(C, ZERO), _mm_andnot_si128(C, ONE));
            _mm_store_si128(&vt->E[i][j], tE);

		}
	}
	
}


