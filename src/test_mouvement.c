#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrdef.h"
#include "nrutil.h"
#include "util.h"
#include "img.h"
#include "mouvement.h"
#include "test_mouvement.h"

void test_corps_boucle1_SigmaDelta_step1() {

	int pass = 1;
	uint8 test;
	uint8 t_1M = 132; // t_1->M[i][j]
	uint8 tI = 86; // t->I[i][j]
	uint8 tM = 131; // t->M[i][j]

	if (t_1M > tI)
		test = t_1M - 1;
	if (tM != test) pass = 0;


	t_1M = 86; // t_1->M[i][j]
	tI = 132; // t->I[i][j]
	tM = 87; // t->M[i][j]

	if (t_1M < tI)
		test = t_1M + 1;
	if (tM != test) pass = 0;

	t_1M = 86; // t_1->M[i][j]
	tI = 86; // t->I[i][j]
	tM = 86; // t->M[i][j]

	if (t_1M == tI)
		test = t_1M;
	if (tM != test) pass = 0;

	if (pass) PASS();
	else FAIL();

}

void test_corps_boucle2_SigmaDelta_Step1() {

	int pass = 1;

	uint8 test;
	uint8 tM = 86;
	uint8 tI = 56;
	uint8 t0 = 30;

	test = abs(tM - tI);
	if (t0 != test) pass = 0;

	test = abs(tI - tM);
	if (t0 != test) pass = 0;

	if (pass) PASS();
	else FAIL();

}

void test_corps_boucle3_SigmaDelta_Step1() {

	int pass = 1;

	uint8 test;
	uint8 t_1V = 132;
	uint8 t0 = 86;
	uint8 tV 133;

	/* N = 2 dans mouvement.h */
	if (t_1V < (N * ))

}

void all_test() {

	test_corps_boucle1_SigmaDelta_step1();
	test_corps_boucle2_SigmaDelta_Step1();

}