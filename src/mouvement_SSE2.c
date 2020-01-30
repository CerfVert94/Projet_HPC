/* ------------------------ */
/* --- mouvement_SSSE3.c --- */
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


/*---------------------------------------------------------------------------------------------*/
void copy_vui8matrix_vui8matrix(vuint8** X, long nrl, long nrh, long vmin, long vmax, vuint8** Y) {
/*---------------------------------------------------------------------------------------------*/
	long i;
	int j;
	for (i = nrl; i <= nrh; i++)
		for (j = vmin; j <= vmax; j++ )
			Y[i][j] = X[i][j];
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SSE(vuint8** I, vuint8** M, vuint8** V, long nrl, long nrh, int v0, int v1) {
/*---------------------------------------------------------------------------------------------*/
	//copy_ui8matrix_ui8matrix(t0->I, t0->nrl, t0->nrh, t0->ncl, t0->nch, t0->M);
	copy_vui8matrix_vui8matrix(I, nrl, nrh, v0, v1, M);
	long i, j;
	for (i = 0; i <= nrh; i++)
		for (j = 0; j <= v1; j++)
			V[i][j] = init_vuint8(1);
}

/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step1_SSE(vuint8** I, vuint8** M_1, vuint8** M, long nrl, long nrh, int v0, int v1) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vI, vM_1, vM, TMP;
	vuint8 C1, C2;

	vuint8 CMP = init_vuint8(127);
	vuint8 ONE = init_vuint8(1);
	vuint8 ZERO = init_vuint8(0);
	vsint8 M_ONE = init_vsint8(-1);

	for(long i = nrl+BORD; i <= nrh-BORD; i++) {
		for(long j = v0+1; j <= v1; j++) {
			vI   = _mm_load_si128((vuint8*) &I[i][j]);
			vM_1 = _mm_load_si128((vuint8*) &M[i][j]);

			vec_cmplt(vM_1, vI, C1, CMP);
			// vI = _mm_sub_epi8(vI, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// C1  = _mm_cmplt_epi8(vM_1, vI);
			// vI = _mm_add_epi8(vI, CMP);
			// vM_1 = _mm_add_epi8(vM_1, CMP);

			C2  = _mm_cmpeq_epi8(vI, vM_1);

		    TMP = _mm_or_si128(_mm_and_si128(C1, ONE), _mm_andnot_si128(C1, M_ONE));
		    TMP = _mm_or_si128(_mm_and_si128(C2, ZERO), _mm_andnot_si128(C2, TMP));
		    vM  = _mm_add_epi8(vM_1, TMP);

		    _mm_store_si128((vuint8*) &M[i][j], vM);

		}
	}
	
}
	

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step2_SSE(vuint8** O, vuint8** M, vuint8** I, long nrl, long nrh, int v0, int v1) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vI, vM, vO;

	vuint8 A, B, C;

	for(long i = nrl+BORD; i <= nrh-BORD; i++) {
		for(long j = v0+1; j <= v1; j++) {
			vI = _mm_load_si128((vuint8*) &I[i][j]);
			vM = _mm_load_si128((vuint8*) &M[i][j]);
			
			vO = vec_subabs(vI, vM);

			_mm_store_si128((vuint8*) &O[i][j], vO);
		}
	}
	
}

/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step3_SSE(vuint8** V, vuint8** V_1, vuint8** O, long nrl, long nrh, int v0, int v1) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 t_1V, vO, NvO, vV, TMP;
	vuint8 C1, C2;

	vuint8 CMP = init_vuint8(127);
	vuint8 ONE = init_vuint8(1);
	vuint8 ZERO = init_vuint8(0);
	vsint8 M_ONE = init_vsint8(-1);

	vuint8 vVmax = init_vuint8(Vmax);
	vuint8 vVmin = init_vuint8(Vmin);

	for(long i = nrl+BORD; i <= nrh-BORD; i++) {
		for(long j = v0+1; j <= v1; j++) {
			t_1V = _mm_load_si128((vuint8*) &V[i][j]);
			vO 	 = _mm_load_si128((vuint8*) &O[i][j]);
			NvO = vO;
			for (int k = 0; k < N; k++) _mm_add_epi8(NvO, vO);

			vec_cmplt(t_1V, NvO, C1, CMP);
			// NvO = _mm_sub_epi8(NvO, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// C1  = _mm_cmplt_epi8(vM_1, NvO);
			// NvO = _mm_add_epi8(NvO, CMP);
			// vM_1 = _mm_add_epi8(vM_1, CMP);

			C2  = _mm_cmpeq_epi8(NvO, t_1V);

		    TMP = _mm_or_si128(_mm_and_si128(C1, ONE), _mm_andnot_si128(C1, M_ONE));
		    TMP = _mm_or_si128(_mm_and_si128(C2, ZERO), _mm_andnot_si128(C2, TMP));
		    vV  = _mm_add_epi8(t_1V, TMP);
		    vV  = _mm_max_epu8(_mm_min_epu8(vV, vVmax), vVmin);

		    _mm_store_si128((vuint8*) &V[i][j], vV);

		}
	}
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step4_SSE(vuint8** O, vuint8** V, vuint8** E, long nrl, long nrh, int v0, int v1) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vO, vV, tE;

	vuint8 ZERO = init_vuint8(0);
	vuint8 ONE = init_vuint8(1);

	vuint8 C, CMP = init_vuint8(127);

	for(long i = nrl+BORD; i <= nrh-BORD; i++) {
		for(long j = v0+1; j <= v1; j++) {
			vO = _mm_load_si128((vuint8*) &O[i][j]);
			vV = _mm_load_si128((vuint8*) &V[i][j]);
			
			vec_cmplt(vO, vV, C, CMP);
			// vO = _mm_sub_epi8(vO, CMP);
			// vV = _mm_sub_epi8(vV, CMP);
			// C  = _mm_cmplt_epi8(vO, vV);
			// vO = _mm_add_epi8(vO, CMP);
			// vV = _mm_add_epi8(vV, CMP);

    		tE = _mm_or_si128(_mm_and_si128(C, ZERO), _mm_andnot_si128(C, ONE));
            _mm_store_si128(&E[i][j], tE);
		}
	}
	
}




