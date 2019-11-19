#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrdef.h"
#include "nrutil.h"
#include "util.h"
#include "img.h"
#include "mouvement.h"
#include "test_mouvement.h"

const char *nom_func;

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
	test = test_corps_SigmaDelta_step1(125, 130);
	UNIT_TEST(test == 126, nom_func, "_msblt");
	/* t_1M > tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step1(130, 127);
	UNIT_TEST(test == 129, nom_func, "_msbgt");
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
	test = test_corps_SigmaDelta_step2(125, 130);
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



