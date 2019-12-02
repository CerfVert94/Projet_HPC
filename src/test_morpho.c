#include <nrdef.h>
#include <nrutil.h>
#include <img.h>
#include <morpho.h>
#include <stdio.h>
#include <stdlib.h>
#include <test_morpho.h>
#include <util.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>
#include <time.h>

const char * nom_func;
#define DEFAULT_PARAMS uint8 **rect, long row, long col, p_struct_elem_dim s

unsigned long long morpho_cpu_cycles(morpho_func_t morpho, uint8 **ppInput, p_struct_elem_dim s, long nrl, long nrh, long ncl, long nch, uint8 **ppOutput)
{
    unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	morpho(ppInput, nrl, nrh, ncl, nch, s, ppOutput);
	end = __rdtsc();
	return end - begin;
}

#define PERM_TO_ROW(perm, irow, ncol) ((perm >> irow * ncol) & 0x1FF) 
#define  ROW_TO_COL(row, icol)        ((row >> (icol)) & 0x1)
void test_morpho()
{
	const long TEST_SE_NROW     = 5L;
	const long TEST_SE_NCOL     = 5L;
	const long TEST_SE_ORIGIN_X = 2L;
	const long TEST_SE_ORIGIN_Y = 2L;
	const long nrl =  TEST_SE_ORIGIN_Y - (TEST_SE_NROW - 1);
	const long nrh = -TEST_SE_ORIGIN_X + (TEST_SE_NROW - 1);
	const long ncl =  TEST_SE_ORIGIN_X - (TEST_SE_NCOL - 1);
	const long nch = -TEST_SE_ORIGIN_X + (TEST_SE_NCOL - 1);

	p_struct_elem_dim s = compute_struct_elem_dim(TEST_SE_ORIGIN_X, TEST_SE_ORIGIN_Y, 
	 											  TEST_SE_NROW    , TEST_SE_NCOL    );
											   
	puts("Checking the structuring element used for unit tests.");
	assert(s->nrow == 5 /* Structuring Element : Height  = 5 */);
	assert(s->ncol == 5 /* Structuring Element : Width   = 5 */);
	assert(s->oriy == 2 /* Structuring Element : originY = 2 */);
	assert(s->orix == 2 /* Structuring Element : originX = 2 */);
	assert(nrl == -2    /* Test input rect : nrl == -2*/);
	assert(nrh ==  2    /* Test input rect : nrh ==  2*/);
	assert(ncl == -2    /* Test input rect : ncl == -2*/);
	assert(nch ==  2    /* Test input rect : nch ==  2*/);
	puts("The structuring element has the correct dimension for tests.");
	// uint8 rect[TEST_SE_NROW][TEST_SE_NCOL];

	uint8 **rect = ui8matrix(nrl, nrh, ncl, nch);
	
	// A binary 5x5 rectangle has 2^25-1 combinations
	// [0 ~ 1FFFFF]
	uint32 min = 0, max = 0x1FFFFFF; // 2^25
	uint32 perm = 0, row = 0;
	const uint32 print_cnt = 0xFFFFF;

	// Loop 1: Get permutations :
	for (perm = 0; perm < max + 1; perm++) {
		// Loop 2 ~ 3 : Tranform the binary permutation to 5x5 rectangle.
		for (int i = nrl; i < nrh + 1; i++) {
			// Get the i-th row from the permutation.
			row = PERM_TO_ROW(perm, (i - nrl), TEST_SE_NCOL); 
			for (int j = ncl; j < nch + 1; j++) {
				// Get the j-th column from the row.
				rect[i][j] = ROW_TO_COL(row, (j - ncl));
			}
		}

		if (perm % print_cnt == 0) 
			printf("%8dth test finished (%08d / %8d).\n", perm, perm, max);
			
		if (perm == max)	assert( Erosion_Test_5x5_Rect(rect, s) == 1 /*perm == 0x1FFFFFF */);
		else				assert( Erosion_Test_5x5_Rect(rect, s) == 0 /*perm <  0x1FFFFFF */);
		if (perm == min)	assert(Dilation_Test_5x5_Rect(rect, s) == 0 /*perm == 0x0000000 */);
		else 				assert(Dilation_Test_5x5_Rect(rect, s) == 1 /*perm >  0x0000000 */);
	}
	puts("Morpho : Passed all tests.");
}
uint8 Erosion_Test_5x5_Rect(uint8 **ppRect, p_struct_elem_dim s)
{
	uint8     output = -1;
	uint8  * pOutput = &output;
	uint8 **ppOutput = &pOutput;
	ui8matrix_erosion_naive(ppRect, 0, 0, 0, 0, s, ppOutput);
	return output;
}

uint8 Dilation_Test_5x5_Rect(uint8 **ppRect, p_struct_elem_dim s)
{
	uint8     output = -1;
	uint8  * pOutput = &output;
	uint8 **ppOutput = &pOutput;
	ui8matrix_dilation_naive(ppRect, 0, 0, 0, 0, s, ppOutput);
	return output;
}