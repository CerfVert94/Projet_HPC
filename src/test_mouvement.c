#include <stdlib.h>
#include <stdio.h>
#include "nrdef.h"
#include "nrutil.h"
#include "util.h"
#include "img.h"
#include "mouvement.h"
#include "test_mouvement.h"

void test_corps_boucle_SigmaDelta_step0() {

	uint8 a; // t0->V[i][j]
	a = 1;
	if (a == 1) PASS();
	else FAIL();

}

void test_corps_boucle1_SigmaDelta_step1() {

	int pass = 1;
	uint8 t_1M = 300; // t_1->M[i][j]
	uint8 tI = 86; // t->I[i][j]
	uint8 tM = 56; // t->M[i][j]
	printf("%d %d\n%d %d\n", t_1M, tI, t_1M-128, tI-128);
	if ((t_1M-128) > (tI-128)) {
		tM = t_1M - 1;	
	}
	if (tM != t_1M - 1) pass = 0;


	t_1M = 86; // t_1->M[i][j]
	tI = 132; // t->I[i][j]
	tM = 56; // t->M[i][j]

	if (t_1M < tI) {
		tM = t_1M + 1;
		if (tM != t_1M + 1) pass = 0;
	}

	t_1M = 86; // t_1->M[i][j]
	tI = 86; // t->I[i][j]
	tM = 56; // t->M[i][j]

	if (t_1M == tI) {
		tM = t_1M;
		if (tM != t_1M) pass = 0;
	}

	if (pass) PASS();
	else FAIL();

}

void test_corps_boucle2_SigmaDelta_Step1() {

	int pass = 1;

}


void all_test() {

	test_corps_boucle_SigmaDelta_step0();
	test_corps_boucle1_SigmaDelta_step1();

}