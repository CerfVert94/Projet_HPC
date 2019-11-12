#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrdef.h"
#include "nrutil.h"
#include "util.h"
#include "img.h"
#include "mouvement.h"
#include "test_mouvement.h"

uint8 test_corps_SigmaDelta_step1(uint8 t_1M, uint8 tI) {

	int pass = 1;

	if (t_1M < tI)
		return t_1M + 1;
	else if (t_1M > tI)
		return t_1M - 1;
	else
		return t_1M;

}

void test_corps_SigmaDelta_step2() {

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

void test_corps_SigmaDelta_step3() {

	int pass = 1;

	uint8 test;
	uint8 t_1V = 132;
	uint8 t0 = 86;
	uint8 tV = 133;

	/* N = 2 dans mouvement.h */
	// if (t_1V < (N * )) 

}

void all_test() {

	uint8 test;

	test = test_corps_SigmaDelta_step1(132, 86);
	if (test == 131) PASS();
	else {FAIL();}

	test = test_corps_SigmaDelta_step1(86, 132);
	if (test == 87) PASS();
	else {FAIL();}

	test = test_corps_SigmaDelta_step1(86, 86);
	if (test == 86) PASS();
	else {FAIL();}

}



