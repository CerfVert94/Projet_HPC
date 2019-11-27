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

/*-----------------------------------------*/
void SigmaDelta_step1_SSE(p_vimage vt, p_vimage vt_1) {
/*-----------------------------------------*/

	vuint8 t, t_1;
	vuint8 vthresh = init_vuint8(THRESHOLD);
	vuint8 vdiff;

	vuint8 A, B, C;

	for(int i = vt->nrl+BORD; i < vt->nrh-BORD; i++) {
		for(int j = vt->v0+1; j < vt->v1-1; j++) {
			t   = _mm_load_si128((vuint8*) &vt->I[i][j]);
			t_1 = _mm_load_si128((vuint8*) &vt_1->I[i][j]);



		}
	}
	
}
	




