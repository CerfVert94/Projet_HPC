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

			vec_cmplsb(tI, t_1I, c, S);
			


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


