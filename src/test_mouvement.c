#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
#include "util.h"
#include "img.h"
#include "mouvement.h"
#include "test_mouvement.h"
#include <stdbool.h>

const char *nom_func;

/*---------------------------------------------------*/
void test_implementation_SigmaDelta_step0(uint8** M, uint8** I, uint8** V, uint8 _Vmin, long nrl, long nrh, long ncl, long nch) {
/*---------------------------------------------------*/
	uint8** V0 = filled_ui8matrix(nrl, nrh, ncl, nch, _Vmin);
	int M_correctly_initialized = !memcmp_ui8matrix(M,  I, nrl, nrh, ncl, nch);
	int V_correctly_initialized = !memcmp_ui8matrix(V, V0, nrl, nrh, ncl, nch);
	// All matrices are correctly initialized.
	assert(M_correctly_initialized && V_correctly_initialized);
}

/*---------------------------------------------------*/
void test_implementation_SigmaDelta_step1(struct sd_set *sd) {
/*---------------------------------------------------*/
	
	uint8** M_t0;
	uint8** I_t0;
	uint8** M_t1;
	long nrl, nrh, ncl, nch;
	sd->sd_func();
}

/*---------------------------------------------------*/
bool SD_step1_produces_valid_output(uint8 m_t0, uint8 i_t1, uint8 m_t1) {
/*---------------------------------------------------*/
	// if P then Q => not P or Q
	int m_lt_i_condition_satisfied = !(m_t0 <  i_t1) || (m_t1 == m_t0 + 1);
	int m_gt_i_condition_satisfied = !(m_t0 >  i_t1) || (m_t1 == m_t0 - 1);
	int m_eq_i_condition_satisfied = !(m_t0 == i_t1) || (m_t1 == m_t0);

	return m_lt_i_condition_satisfied &&
	 	   m_gt_i_condition_satisfied &&
	 	   m_eq_i_condition_satisfied;
}

/*---------------------------------------------------*/
bool SD_step2_produces_valid_output(uint8 m_t, uint8 i_t, uint8 o_t) {
/*---------------------------------------------------*/
	return (abs(m_t - i_t) == o_t);
}
/*---------------------------------------------------*/
bool SD_step3_produces_valid_output(uint8 v_t0, uint8 n, uint8 o_t, uint8 v_t1, uint8 _vmin, uint8 _vmax) {
/*---------------------------------------------------*/
	int min_clamped = (0 < v_t0       && v_t0 - 1 < _vmin && v_t1 == _vmin) ||
					  (                  v_t0 - 0 < _vmin && v_t1 == _vmin);
	int max_clamped = (    v_t0 < 255 && v_t0 + 1 > _vmax && v_t1 == _vmax) ||
					  (                  v_t0 + 0 > _vmax && v_t1 == _vmax);
	// if P then Q => not P or Q
	int v_lt_N_o_condition_satisfied = !(v_t0 <  n * o_t) || ((v_t1 == v_t0 + 1)                || max_clamped);
	int v_gt_N_o_condition_satisfied = !(v_t0 >  n * o_t) || ((v_t1 == v_t0 - 1) || min_clamped               );
	int v_eq_N_o_condition_satisfied = !(v_t0 == n * o_t) || ((v_t1 == v_t0)     || min_clamped || max_clamped);
	

	return v_lt_N_o_condition_satisfied && 
		   v_gt_N_o_condition_satisfied &&
	 	   v_eq_N_o_condition_satisfied;
}

/*---------------------------------------------------*/
bool SD_step4_produces_valid_output(uint8 o_t,uint8 v_t, uint8 e) {
/*---------------------------------------------------*/
	int o_lt_v_condition_satisfied  = !(o_t <  v_t) || (e == 0);
	int o_geq_v_condition_satisfied = !(o_t >= v_t) || (e == 1);

	return o_lt_v_condition_satisfied && 
	       o_geq_v_condition_satisfied;
}
/*---------------------------------------------------*/
uint8 test_corps_SigmaDelta_step1(uint8 t_1M, uint8 tI) {
/*---------------------------------------------------*/
	nom_func = __func__;

	if (t_1M < tI)
		return t_1M + 1;
	else if (t_1M > tI)
		return t_1M - 1;
	else
		return t_1M;

}

/*-------------------------------------------------*/
uint8 test_corps_SigmaDelta_step2(uint8 tM, uint8 tI) {
/*-------------------------------------------------*/
	nom_func = __func__;

	return abs(tM - tI);

}

/*---------------------------------------------------*/
uint8 test_corps_SigmaDelta_step3(uint8 t_1V, uint8 tO) { 
/*---------------------------------------------------*/

	int M = 2;

 	nom_func = __func__;

	uint8 tV;

	if (t_1V < (M * tO))
		tV = t_1V + 1;
	else if (t_1V > (M * tO))
		tV = t_1V - 1;
	else
		tV = t_1V;
	return max(min(tV, Vmax), Vmin);

}

/*---------------------------------------------------*/
uint8 test_corps_SigmaDelta_step4(uint8 tO, uint8 tV) {
/*---------------------------------------------------*/

	nom_func = __func__;

	if ( tO < tV )
		return 0;
	else
		return 1;

}

/*---------------------*/
void all_test_mouvement() { /* UNIT_TEST dans util.h */
/*---------------------*/
	uint8 test;

	printf("\n==== test_corps_SigmaDelta_step1 ====\n");
	/* t_1M < tI */
	test = test_corps_SigmaDelta_step1(15, 30);
	UNIT_TEST(test == 16, nom_func, "_lt");
	/* t_1M > tI */
	test = test_corps_SigmaDelta_step1(30, 15);
	UNIT_TEST(test == 29, nom_func, "_gt");
	/* t_1M = tI */
	test = test_corps_SigmaDelta_step1(0, 0);
	/* t_1M > */
	UNIT_TEST(test == 0, nom_func, "_eq");
	/* t_1M < tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step1(127, 128);
	UNIT_TEST(test == 128, nom_func, "_msblt");
	/* t_1M > tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step1(128, 127);
	UNIT_TEST(test == 127, nom_func, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step2 ====\n");
	/* tM < tI */
	test = test_corps_SigmaDelta_step2(15, 30);
	UNIT_TEST(test == 15, nom_func, "_lt");
	/* tM > tI */
	test = test_corps_SigmaDelta_step2(30, 15);
	UNIT_TEST(test == 15, nom_func, "_gt");
	/* tM = tI */
	test = test_corps_SigmaDelta_step2(0, 0);
	UNIT_TEST(test == 0, nom_func, "_eq");
	/* tM < tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step2(127, 128);
	UNIT_TEST(test == 5, nom_func, "_msblt");
	/* tM > tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step2(130, 125);
	UNIT_TEST(test == 5, nom_func, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step3 ====\n");
	/* t_1V < tI */
	test = test_corps_SigmaDelta_step3(14, 30);
	UNIT_TEST(test == 15, nom_func, "_lt");
	/* t_1V > tI */
	test = test_corps_SigmaDelta_step3(30, 14);
	UNIT_TEST(test == 29, nom_func, "_ht");
	/* t_1V = tI */
	test = test_corps_SigmaDelta_step3(0, 0);
	UNIT_TEST(test == 1, nom_func, "_eq");
	/* t_1V < tI compare signe */
	test = test_corps_SigmaDelta_step3(61, 130);
	UNIT_TEST(test == 62, nom_func, "_msblt");
	test = test_corps_SigmaDelta_step3(130, 61);
	UNIT_TEST(test == 129, nom_func, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step4 ====\n");
	/* t_0 < tV */
	test = test_corps_SigmaDelta_step3(14, 30);
	UNIT_TEST(test == 15, nom_func, "_lt");
	/* t_1V > tI */
	test = test_corps_SigmaDelta_step3(30, 14);
	UNIT_TEST(test == 29, nom_func, "_ht");
	/* t_1V = tI */
	test = test_corps_SigmaDelta_step3(0, 0);
	UNIT_TEST(test == 1, nom_func, "_eq");
	/* t_1V < tI compare signe */
	test = test_corps_SigmaDelta_step3(61, 130);
	UNIT_TEST(test == 62, nom_func, "_msblt");
	test = test_corps_SigmaDelta_step3(130, 61);
	UNIT_TEST(test == 129, nom_func, "_msbgt");
	printf("=====================================\n\n");


}



