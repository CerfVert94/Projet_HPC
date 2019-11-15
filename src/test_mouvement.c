#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrdef.h"
#include "nrutil.h"
#include "util.h"
#include "img.h"
#include "mouvement.h"
#include "test_mouvement.h"

char *nom_func;

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
void all_test_mouvement() { /* UTEST dans util.h */
/*---------------------*/
	uint8 test;

	printf("\n==== test_corps_SigmaDelta_step1 ====\n");
	/* t_1M < tI */
	test = test_corps_SigmaDelta_step1(15, 30);
	UTEST(test == 16, "_lt");
	/* t_1M > tI */
	test = test_corps_SigmaDelta_step1(30, 15);
	UTEST(test == 29, "_ht");
	/* t_1M = tI */
	test = test_corps_SigmaDelta_step1(0, 0);
	/* t_1M > */
	UTEST(test == 0, "_eq");
	/* t_1M < tI compare signe */
	test = test_corps_SigmaDelta_step1(125, 130);
	UTEST(test == 126, "_msblt");
	/* t_1M > tI compare signe */
	test = test_corps_SigmaDelta_step1(130, 127);
	UTEST(test == 129, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step2 ====\n");
	/* tM < tI */
	test = test_corps_SigmaDelta_step2(15, 30);
	UTEST(test == 15, "_lt");
	/* tM > tI */
	test = test_corps_SigmaDelta_step2(30, 15);
	UTEST(test == 15, "_gt");
	/* tM = tI */
	test = test_corps_SigmaDelta_step2(0, 0);
	UTEST(test == 0, "_eq");
	/* tM < tI compare signe */
	test = test_corps_SigmaDelta_step2(125, 130);
	UTEST(test == 5, "_msblt");
	/* tM > tI compare signe */
	test = test_corps_SigmaDelta_step2(130, 125);
	UTEST(test == 5, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step3 ====\n");
	/* t_1V < tI */
	test = test_corps_SigmaDelta_step3(14, 30);
	UTEST(test == 15, "_lt");
	/* t_1V > tI */
	test = test_corps_SigmaDelta_step3(30, 14);
	UTEST(test == 29, "_ht");
	/* t_1V = tI */
	test = test_corps_SigmaDelta_step3(0, 0);
	UTEST(test == 1, "_eq");
	/* t_1V < tI compare signe */
	test = test_corps_SigmaDelta_step3(61, 130);
	UTEST(test == 62, "_msblt");
	test = test_corps_SigmaDelta_step3(130, 61);
	UTEST(test == 129, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step4 ====\n");
	/* t_0 < tV */
	test = test_corps_SigmaDelta_step3(14, 30);
	UTEST(test == 15, "_lt");
	/* t_1V > tI */
	test = test_corps_SigmaDelta_step3(30, 14);
	UTEST(test == 29, "_ht");
	/* t_1V = tI */
	test = test_corps_SigmaDelta_step3(0, 0);
	UTEST(test == 1, "_eq");
	/* t_1V < tI compare signe */
	test = test_corps_SigmaDelta_step3(61, 130);
	UTEST(test == 62, "_msblt");
	test = test_corps_SigmaDelta_step3(130, 61);
	UTEST(test == 129, "_msbgt");
	printf("=====================================\n\n");


}



