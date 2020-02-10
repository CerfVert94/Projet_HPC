/* ------------------------ */
/* --- mouvement_SIMD.c --- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <tmmintrin.h>
#include <pmmintrin.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"
#include "myvnrutil.h"

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "img_SIMD.h"
#include "mouvement_SIMD.h"
#include "util.h"



/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SIMD(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	long i, j;


	vuint8 vevLT = _mm_set_epi8(v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min);
	vuint8 vec0 = _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	vuint8 *ptr = &V[nrl][v0];
	
	memcpy(&M[nrl][v0], &I[nrl][v0], (nrh - nrl + 1) * (v1 - v0 + 1) *sizeof(vuint8));
	for (i = nrl; i < nrh + 1; i++)
		for (j = v0; j < v1 + 1; j++)
			V[i][j] = vevLT;
}
/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SIMD_and_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	long i, j;
	vuint8 vevLT = _mm_set_epi8(v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min);
	vuint8 vec0 = _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	vuint8 *ptr = &V[nrl][v0];
	
	memcpy(&M[nrl][v0], &I[nrl][v0], (nrh - nrl + 1) * (v1 - v0 + 1) *sizeof(vuint8));
	for (i = nrl; i < nrh + 1; i++)
		for (j = v0; j < v1 + 1; j++)
			V[i][j] = _mm_and_si128(vevLT, vevLT);
}
/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SIMD_or_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	long i, j;
	vuint8 vevLT = _mm_set_epi8(v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min);
	
	memcpy(&M[nrl][v0], &I[nrl][v0], (nrh - nrl + 1) * (v1 - v0 + 1) *sizeof(vuint8));
	for (i = nrl; i < nrh + 1; i++)
		for (j = v0; j < v1 + 1; j++)
			V[i][j] = _mm_or_si128(vevLT, vevLT);
}
/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SIMD_load_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	long i, j;
	vuint8 vevLT = _mm_set_epi8(v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min);
	vuint8 vec0 = _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	vuint8 *ptr = &V[nrl][v0];
	
	memcpy(&M[nrl][v0], &I[nrl][v0], (nrh - nrl + 1) * (v1 - v0 + 1) *sizeof(vuint8));
	for (i = nrl; i < nrh + 1; i++)
		for (j = v0; j < v1 + 1; j++)
			V[i][j] = _mm_lddqu_si128(&vevLT);
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SIMD_store_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	long i, j;
	vuint8 vevLT = _mm_set_epi8(v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min, v_min);
	vuint8 vec0 = _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	vuint8 *ptr = &V[nrl][v0];
	
	memcpy(&M[nrl][v0], &I[nrl][v0], (nrh - nrl + 1) * (v1 - v0 + 1) *sizeof(vuint8));
	for (i = nrl; i < nrh + 1; i++)
		for (j = v0; j < v1 + 1; j++)
			_mm_store_si128(&V[i][j], vevLT);
}
/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SIMD_memset_load(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	long i, j;
	uint32 one = 1 << 0 | 1 << 8 | 1 << 16 | 1 << 24;
	vuint8 *ptr = &V[nrl][v0];
	size_t len = (nrh - nrl + 1) * (v1 - v0 + 1)  * sizeof(vuint8);
	
	memcpy(&M[nrl][v0], &I[nrl][v0], len);
	memset((uint32*)&V[nrl][v0], one, len);
}





/*-----------------------------------------------------------------------------------------------*/
vuint8 SigmaDelta_step1_SIMD_single_vector(vuint8 vM_1, vuint8 vI, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vRes, vTemp;
	vuint8 vLT, vEQ, C3;

	vuint8 CMP = _mm_set1_epi8(128);
	vuint8 ONE = _mm_set1_epi8(0x1);
	vuint8 vM_1sub, vM_1add;

	vM_1sub = _mm_subs_epu8(vM_1, ONE);
	vM_1add = _mm_adds_epu8(vM_1, ONE);

	vec_cmplt(vM_1, vI, vLT, CMP);
	vEQ  = _mm_cmpeq_epi8(vI, vM_1);
	// [vM_1 < vI => vTemp := vM_1add] or [not (vM_1 < vI) => vTemp := vM_1add]
	vTemp = _mm_or_si128(_mm_and_si128(vLT, vM_1add), _mm_andnot_si128(vLT, vM_1sub));
	// If vM1 == vI, then vM := vM_1. If not vM := vTemp (less than / greater than )
	vRes  = _mm_or_si128(_mm_and_si128(vEQ, vM_1)   , _mm_andnot_si128(vEQ, vTemp));
	return vRes;

	
}
/*---------------------------------------------------------------------------------------------*/
vuint8 SigmaDelta_step2_SIMD_single_vector(vuint8 vM, vuint8 vI, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	return vec_subabs(vI, vM);	
}


/*-----------------------------------------------------------------------------------------------*/
vuint8 SigmaDelta_step3_SIMD_single_vector(vuint8 vV_1, vuint8 vO, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vNO, vTemp, vRes;
	vuint8 vLT, vEQ;

	vuint8 CMP = _mm_set1_epi8(128);
	vuint8 ONE = _mm_set1_epi8(0x1);
	vuint8 vV_1add, vV_1sub;


	vuint8 vuint8_max = _mm_set1_epi8(255);
	vuint8 vVmax = _mm_set1_epi8(v_max);
	vuint8 vVmin = _mm_set1_epi8(v_min);

	vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
	vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
	
	vNO = vO;
	for (int k = 1; k < n_coeff; k++) vNO =_mm_adds_epu8(vNO, vO);

	// V_t0 < n * O_t1
	vec_cmplt(vV_1, vNO, vLT, CMP);
	// V_t0 == n * O_t1
	vEQ  = _mm_cmpeq_epi8(vNO, vV_1);

	vTemp = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_andnot_si128(vLT, vV_1sub));
	vTemp = _mm_or_si128(_mm_and_si128(vEQ, vV_1)   , _mm_andnot_si128(vEQ, vTemp));
	vRes  = _mm_max_epu8(_mm_min_epu8(vTemp, vVmax), vVmin);
	return vRes;
}
/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step3_SIMD_ver2(vuint8** V_1, vuint8** O, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vV_1, vO, vNO, vV, vTemp;
	vuint8 vLT, vEQ, C3;

	vuint8 CMP = _mm_set1_epi8(128);
	vuint8 ONE = _mm_set1_epi8(1);
	vuint8 vV_1add, vV_1sub;
	vuint8 vVmax = _mm_set1_epi8(v_max);
	vuint8 vVmin = _mm_set1_epi8(v_min);
	

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vV_1 = _mm_load_si128((vuint8*) &V_1[i][j]); 
			vO 	 = _mm_load_si128((vuint8*) &O[i][j]);   

			vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
			vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
			
			vNO = vO;
			for (int k = 1; k < n_coeff; k++) vNO =_mm_adds_epu8(vNO, vO);

			// V_t0 < n * O_t1
			vec_cmplt(vV_1, vNO, vLT, CMP);
			// V_t0 == n * O_t1
			vEQ  = _mm_cmpeq_epi8(vNO, vV_1);


			vTemp = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_andnot_si128(vLT, vV_1sub));
			vTemp = _mm_or_si128(_mm_and_si128(vEQ, vV_1)   , _mm_andnot_si128(vEQ, vTemp));
			vV  = _mm_max_epu8(_mm_min_epu8(vTemp, vVmax), vVmin);

		    _mm_store_si128((vuint8*) &V[i][j], vV);

		}
	}
}


/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step1_SIMD(vuint8** M_1, vuint8** I, vuint8** M, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vI, vM_1, vM, TMP, A, B;
	vuint8 vLT, vEQ, C3;

	vuint8 CMP = init_vuint8(128);
	vuint8 ONE = init_vuint8(1);
	vuint8 vM_1sub, vM_1add;

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vI   = _mm_load_si128((vuint8*) &I[i][j]);
			vM_1 = _mm_load_si128((vuint8*) &M_1[i][j]);
			vM_1sub = _mm_subs_epu8(vM_1, ONE);
			vM_1add = _mm_adds_epu8(vM_1, ONE);

			vec_cmplt(vM_1, vI, vLT, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// vI = _mm_sub_epi8(vI, CMP);
			// vLT = _mm_cmplt_epi8(vM_1, vI);
			// vM_1 = _mm_add_epi8(vM_1, CMP);
			// vI = _mm_add_epi8(vI, CMP);
			vEQ  = _mm_cmpeq_epi8(vI, vM_1);
		    TMP = _mm_or_si128(_mm_and_si128(vLT, vM_1add), _mm_andnot_si128(vLT, vM_1sub));
		    vM = _mm_or_si128(_mm_and_si128(vEQ, vM_1), _mm_andnot_si128(vEQ, TMP));
		    _mm_store_si128((vuint8*) &M[i][j], vM);

		}
	}
	
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step2_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vI, vM, vO;
	vuint8 vMI_sub, vShift;

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vI = _mm_load_si128((vuint8*) &I[i][j]);
			vM = _mm_load_si128((vuint8*) &M[i][j]);
			
			vO = vec_subabs(vI, vM);
			
			_mm_store_si128((vuint8*) &O[i][j], vO);
		}
	}
	
}
/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step2_ExLU_O3_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max){
/*---------------------------------------------------------------------------------------------*/
	vuint8 vI, vM, vO;
	vuint8 *vI_row0, *vI_row1, *vI_row2;
	vuint8 *vM_row0, *vM_row1, *vM_row2;
	vuint8 *vO_row0, *vO_row1, *vO_row2;

	vuint8 vMI_sub, vShift;
	long order = 3;
	long r = (nrh - nrl + 1) % order;
	
	for(long i = nrl; i <= nrh - r; i += order) {
		vI_row0 = I[i+0]; vI_row1 = I[i+1]; vI_row2 = I[i+2];
		vM_row0 = M[i+0]; vM_row1 = M[i+1]; vM_row2 = M[i+2];
		vO_row0 = O[i+0]; vO_row1 = O[i+1]; vO_row2 = O[i+2];
		for(long j = v0; j <= v1; j++) {
			vI = _mm_load_si128((vuint8*) &vI_row0[j]);
			vM = _mm_load_si128((vuint8*) &vM_row0[j]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row0[j], vO);

			vI = _mm_load_si128((vuint8*) &vI_row1[j]);
			vM = _mm_load_si128((vuint8*) &vM_row1[j]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row1[j], vO);

			vI = _mm_load_si128((vuint8*) &vI_row2[j]);
			vM = _mm_load_si128((vuint8*) &vM_row2[j]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row2[j], vO);
		}
	}
	switch(r){
		case 2:
			vI_row0 = I[nrh - 1]; vI_row1 = I[nrh]; 
			vM_row0 = M[nrh - 1]; vM_row1 = M[nrh]; 
			vO_row0 = O[nrh - 1]; vO_row1 = O[nrh]; 
			for(long j = v0; j <= v1; j++) {
				vI = _mm_load_si128((vuint8*) &vI_row0[j]);
				vM = _mm_load_si128((vuint8*) &vM_row0[j]);
				vO = vec_subabs(vI, vM);
				_mm_store_si128((vuint8*) &vO_row0[j], vO);

				vI = _mm_load_si128((vuint8*) &vI_row1[j]);
				vM = _mm_load_si128((vuint8*) &vM_row1[j]);
				vO = vec_subabs(vI, vM);
				_mm_store_si128((vuint8*) &vO_row1[j], vO);
			}
			break;
		case 1:
			vI_row1 = I[nrh]; 
			vM_row1 = M[nrh]; 
			vO_row1 = O[nrh]; 
			for(long j = v0; j <= v1; j++) {
				vI = _mm_load_si128((vuint8*) &vI_row1[j]);
				vM = _mm_load_si128((vuint8*) &vM_row1[j]);
				vO = vec_subabs(vI, vM);
				_mm_store_si128((vuint8*) &vO_row1[j], vO);
			}
			break;
	}
	
}
/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step2_InLU_O3_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vI, vM, vO;
	vuint8  *vI_row, *vM_row, *vO_row;
	vuint8 vMI_sub, vShift;
	long order = 3;
	long r = (v1 - v0 + 1) % order;

	for(long i = nrl; i <= nrh; i++) {
		vI_row = I[i];
		vM_row = M[i];
		vO_row = O[i];
		for(long j = v0; j <= v1 - r; j+= order) {
			vI = _mm_load_si128((vuint8*) &vI_row[j + 0]);
			vM = _mm_load_si128((vuint8*) &vM_row[j + 0]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row[j+0], vO);

			vI = _mm_load_si128((vuint8*) &vI_row[j + 1]);
			vM = _mm_load_si128((vuint8*) &vM_row[j + 1]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row[j+1], vO);
			
			vI = _mm_load_si128((vuint8*) &vI_row[j + 2]);
			vM = _mm_load_si128((vuint8*) &vM_row[j + 2]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row[j]+2, vO );
		}
	}
	switch(r) {
		case 2:

		for(long i = nrl; i <= nrh; i++) {
			
			vI = _mm_load_si128((vuint8*) &I[i][v1 - 1]);
			vM = _mm_load_si128((vuint8*) &M[i][v1 - 1]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &O[i][v1 -1], vO);

			vI = _mm_load_si128((vuint8*) &I[i][v1 + 0]);
			vM = _mm_load_si128((vuint8*) &M[i][v1 + 0]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &O[i][v1 - 0], vO);
		}
		break;
		case 1:

		for(long i = nrl; i <= nrh; i++) {
			vI = _mm_load_si128((vuint8*) &I[i][v1 + 0]);
			vM = _mm_load_si128((vuint8*) &M[i][v1 + 0]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &O[i][v1 - 0], vO);
		}
		break;
		case 0:
		break;
	}
	
}




/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step3_SIMD(vuint8** V_1, vuint8** O, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vV_1, vO, NvO, vV, TMP;
	vuint8 vLT, vEQ, C3;

	vuint8 CMP = init_vuint8(128);
	vuint8 ONE = init_vuint8(1);
	vuint8 vV_1add, vV_1sub;


	vuint8 vuint8_max = init_vuint8(255);
	vuint8 vVmax = init_vuint8(v_max);
	vuint8 vVmin = init_vuint8(v_min);
	// vuint8 vVmin = init_vuint8(v_min);
	// printf("vmin:%u\n",v_min);
	// printf("vmax:%u\n",v_max);

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vV_1 = _mm_load_si128((vuint8*) &V_1[i][j]); 
			vO 	 = _mm_load_si128((vuint8*) &O[i][j]);   

			vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
			vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
			// display_vuint8(vV_1, " %03u", "vV_1\t\t"); 			printf("\n");
			// display_vuint8(vV_1sub, " %03u", "vV_1sub\t\t");	printf("\n");
			// display_vuint8(vV_1add, " %03u", "vV_1add\t\t");	printf("\n");
			NvO = vO;
			for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
			// display_vuint8(NvO, " %03u", "NvO\t\t");	printf("\n");

			// V_t0 < n * O_t1
			vec_cmplt(vV_1, NvO, vLT, CMP);

			// V_t0 == n * O_t1
			vEQ  = _mm_cmpeq_epi8(NvO, vV_1);

			// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
			C3 = _mm_andnot_si128(vLT, _mm_andnot_si128(vEQ, vuint8_max));
			// display_vuint8(vLT, "%4u", "less than\t");	printf("\n");
			// display_vuint8(vEQ, "%4u", "equal to\t");	printf("\n");
			// display_vuint8(C3, "%4u", "greater than\t");printf("\n");
			// NvO = _mm_sub_epi8(NvO, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// vLT  = _mm_cmplt_epi8(vM_1, NvO);
			// NvO = _mm_add_epi8(NvO, CMP);
			// vM_1 = _mm_add_epi8(vM_1, CMP);


			// If V_t0 < n * O_t1 then vV := vV_1add 
			// OR
			// If V_t0 > n * O_t1 then vV := vV_1sub
		    TMP = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_and_si128(C3, vV_1sub));
			// display_vuint8(TMP, "%4u", "TMP\t\t");printf("\n");
		    TMP = _mm_or_si128(TMP, _mm_and_si128(vEQ, vV_1));
			// display_vuint8(TMP, "%4u", "TMP\t\t");printf("\n");
		    vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);

		    _mm_store_si128((vuint8*) &V[i][j], vV);

		}
	}
}



/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step3_InLU_O3_SIMD(vuint8** V_1, vuint8** O, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/	
	long order = 3;
	long r = (v1 - v0 + 1) % order;vuint8 vV_1, vO, NvO, vV, TMP;
	vuint8 vLT, vEQ, C3;

	vuint8 CMP = init_vuint8(128);
	vuint8 ONE = init_vuint8(1);
	vuint8 vV_1add, vV_1sub;


	vuint8 vuint8_max = init_vuint8(255);
	vuint8 vVmax = init_vuint8(v_max);
	vuint8 vVmin = init_vuint8(v_min);
	

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1 - r; j += order) {
			vV_1 = _mm_load_si128((vuint8*) &V_1[i][j]); 
			vO 	 = _mm_load_si128((vuint8*) &O[i][j]);   

			vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
			vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
			
			NvO = vO;
			for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
			
			// V_t0 < n * O_t1
			vec_cmplt(vV_1, NvO, vLT, CMP);
			// V_t0 == n * O_t1
			vEQ  = _mm_cmpeq_epi8(NvO, vV_1);
			// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
			C3 = _mm_andnot_si128(vLT, _mm_andnot_si128(vEQ, vuint8_max));
			
			// If V_t0 > n * O_t1 then vV := vV_1sub
		    TMP = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_and_si128(C3, vV_1sub));
			TMP = _mm_or_si128(TMP, _mm_and_si128(vEQ, vV_1));
			vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);
		    _mm_store_si128((vuint8*) &V[i][j], vV);

			vV_1 = _mm_load_si128((vuint8*) &V_1[i][j + 1]); 
			vO 	 = _mm_load_si128((vuint8*) &O[i][j + 1]);   

			vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
			vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
			
			NvO = vO;
			for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
			
			// V_t0 < n * O_t1
			vec_cmplt(vV_1, NvO, vLT, CMP);
			// V_t0 == n * O_t1
			vEQ  = _mm_cmpeq_epi8(NvO, vV_1);
			// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
			C3 = _mm_andnot_si128(vLT, _mm_andnot_si128(vEQ, vuint8_max));
			
			// If V_t0 > n * O_t1 then vV := vV_1sub
		    TMP = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_and_si128(C3, vV_1sub));
			TMP = _mm_or_si128(TMP, _mm_and_si128(vEQ, vV_1));
			vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);

		    _mm_store_si128((vuint8*) &V[i][j + 1], vV);

			vV_1 = _mm_load_si128((vuint8*) &V_1[i][j+2]); 
			vO 	 = _mm_load_si128((vuint8*) &O[i][j+2]);   

			vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
			vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
			
			NvO = vO;
			for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
			
			// V_t0 < n * O_t1
			vec_cmplt(vV_1, NvO, vLT, CMP);
			// V_t0 == n * O_t1
			vEQ  = _mm_cmpeq_epi8(NvO, vV_1);
			// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
			C3 = _mm_andnot_si128(vLT, _mm_andnot_si128(vEQ, vuint8_max));
			
			// If V_t0 > n * O_t1 then vV := vV_1sub
		    TMP = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_and_si128(C3, vV_1sub));
			TMP = _mm_or_si128(TMP, _mm_and_si128(vEQ, vV_1));
			vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);

		    _mm_store_si128((vuint8*) &V[i][j+2], vV);

		}
	}

	switch(r){
		case 2: 
			for(long i = nrl; i <= nrh; i++) {
				vV_1 = _mm_load_si128((vuint8*) &V_1[i][v1 - 1]); 
				vO 	 = _mm_load_si128((vuint8*) &O[i][v1 - 1]);   

				vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
				vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
				
				NvO = vO;
				for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
				
				// V_t0 < n * O_t1
				vec_cmplt(vV_1, NvO, vLT, CMP);
				// V_t0 == n * O_t1
				vEQ  = _mm_cmpeq_epi8(NvO, vV_1);
				// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
				C3 = _mm_andnot_si128(vLT, _mm_andnot_si128(vEQ, vuint8_max));
				
				// If V_t0 > n * O_t1 then vV := vV_1sub
				TMP = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_and_si128(C3, vV_1sub));
				TMP = _mm_or_si128(TMP, _mm_and_si128(vEQ, vV_1));
				vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);
				_mm_store_si128((vuint8*) &V[i][v1 - 1], vV);

				vV_1 = _mm_load_si128((vuint8*) &V_1[i][v1]); 
				vO 	 = _mm_load_si128((vuint8*) &O[i][v1]);   

				vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
				vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
				
				NvO = vO;
				for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
				
				// V_t0 < n * O_t1
				vec_cmplt(vV_1, NvO, vLT, CMP);
				// V_t0 == n * O_t1
				vEQ  = _mm_cmpeq_epi8(NvO, vV_1);
				// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
				C3 = _mm_andnot_si128(vLT, _mm_andnot_si128(vEQ, vuint8_max));
				
				// If V_t0 > n * O_t1 then vV := vV_1sub
				TMP = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_and_si128(C3, vV_1sub));
				TMP = _mm_or_si128(TMP, _mm_and_si128(vEQ, vV_1));
				vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);

				_mm_store_si128((vuint8*) &V[i][v1], vV);
			}
		case 1: 
			for(long i = nrl; i <= nrh; i++) {

				vV_1 = _mm_load_si128((vuint8*) &V_1[i][v1]); 
				vO 	 = _mm_load_si128((vuint8*) &O[i][v1]);   

				vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
				vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
				
				NvO = vO;
				for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
				
				// V_t0 < n * O_t1
				vec_cmplt(vV_1, NvO, vLT, CMP);
				// V_t0 == n * O_t1
				vEQ  = _mm_cmpeq_epi8(NvO, vV_1);
				// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
				C3 = _mm_andnot_si128(vLT, _mm_andnot_si128(vEQ, vuint8_max));
				
				// If V_t0 > n * O_t1 then vV := vV_1sub
				TMP = _mm_or_si128(_mm_and_si128(vLT, vV_1add), _mm_and_si128(C3, vV_1sub));
				TMP = _mm_or_si128(TMP, _mm_and_si128(vEQ, vV_1));
				vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);

				_mm_store_si128((vuint8*) &V[i][v1], vV);
			}
			break;
		

	}
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step4_SIMD(vuint8** O, vuint8** V, vuint8** E, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vO, vV, vE;

	vuint8 ZERO = init_vuint8(0);
	vuint8 ONE = init_vuint8(1);

	vuint8 C, CMP = init_vuint8(128);

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vO = _mm_load_si128((vuint8*) &O[i][j]);
			vV = _mm_load_si128((vuint8*) &V[i][j]);
			
			vec_cmplt(vO, vV, C, CMP);
			// vO = _mm_sub_epi8(vO, CMP);
			// vV = _mm_sub_epi8(vV, CMP);
			// C  = _mm_cmplt_epi8(vO, vV);
			// vO = _mm_add_epi8(vO, CMP);
			// vV = _mm_add_epi8(vV, CMP);

    		vE = _mm_or_si128(_mm_and_si128(C, ZERO), _mm_andnot_si128(C, ONE));
            _mm_store_si128(&E[i][j], vE);
		}
	}
	
}

/*-----------------------------------------------------------------------------------------------*/
vuint8 SigmaDelta_step4_SIMD_single_vector(vuint8 vO, vuint8 vV, uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/
	vuint8 ZERO = _mm_setzero_si128();
	vuint8 ONE = _mm_set1_epi8(1);

	vuint8 vLT, CMP = _mm_set1_epi8(128);
	vuint8 vRes;

	vec_cmplt(vO, vV, vLT, CMP);

	vRes = _mm_or_si128(_mm_and_si128(vLT, ZERO), _mm_andnot_si128(vLT, ONE));
	return vRes;
}
	
void SigmaDelta_SIMD(p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	SigmaDelta_step1_SIMD(t0->M, t1->I, t1->M, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	SigmaDelta_step2_InLU_O3_SIMD(t1->M, t1->I, t1->O, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	// SigmaDelta_step2_SIMD(t1->M, t1->I, t1->O, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	// SigmaDelta_step3_SIMD(t0->V, t1->O, t1->V, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	SigmaDelta_step3_InLU_O3_SIMD(t0->V, t1->O, t1->V, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	SigmaDelta_step4_SIMD(t1->O, t1->V, t1->E, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
}

void SigmaDelta_SIMD_FL(p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{

	int nrl, nrh, v0, v1;

	vuint8 **M_1 = t0->M, **I_1 = t0->I, **O_1 = t0->O, **V_1 = t0->V, **E_1 = t0->E;
	vuint8 **M   = t1->M, **I   = t1->I, **O   = t1->O, **V   = t1->V, **E   = t1->E;

	nrl = t1->nrl; v0 = t1->v0;
	nrh = t1->nrh; v1 = t1->v1;	

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			M[i][j] = SigmaDelta_step1_SIMD_single_vector(M_1[i][j], I[i][j], n_coeff, v_min, v_max);
			O[i][j] = SigmaDelta_step2_SIMD_single_vector(M  [i][j], I[i][j], n_coeff, v_min, v_max);
			V[i][j] = SigmaDelta_step3_SIMD_single_vector(V_1[i][j], O[i][j], n_coeff, v_min, v_max);
			E[i][j] = SigmaDelta_step4_SIMD_single_vector(O  [i][j], V[i][j], n_coeff, v_min, v_max);
		}
	}
}
void test_step0_SIMD() {

	// p_vimage vT_1 = create_vimage("../car3/car_3000.pgm");
	// p_vimage vT   = create_vimage("../car3/car_3001.pgm");

	// SigmaDelta_step0_SIMD(vT_1->I, vT_1->M, vT_1->V, vT_1->nrl, vT_1->nrh, vT_1->v0, vT_1->v1);
	// SigmaDelta_step1_SIMD(vT->I, vT_1->M, vT->M, vT->nrl, vT->nrh, vT->v0, vT->v1);
}


