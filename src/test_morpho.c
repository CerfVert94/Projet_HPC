#include <nrdef.h>
#include <nrutil.h>
#include <img.h>
#include <morpho.h>
#include <stdio.h>
#include <stdlib.h>
#include <test_morpho.h>
#include <util.h>
#include <mutil.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>
#include <time.h>
uint8 **ui8matrix_test_input(long size, const uint8 xor_mask, p_struct_elem_dim s)
{
	uint8** ppOutput, cell = 255, first_cell;
	const uint8 CELL = 1;
	long nrl, nrh, ncl, nch;
	nrl = 0 + s->nrl;
	nrh = size + s->nrh - 1;
	ncl = 0 + s->ncl;
	nch = size + s->nch - 1;
	ppOutput = ui8matrix(nrl, nrh, ncl, nch);
	memset(&ppOutput[nrl][ncl], 0, NROW(nrh, nrl) * NCOL(nch, ncl) * sizeof(ppOutput[0][0]));
	first_cell = xor_mask;
	for (long i = 0; i < size; i++) {
		if (i % 3 == 0 && i > 0)
			first_cell ^= xor_mask;    
		cell = first_cell;
		for (long j = 0; j < size; j++) {
			if (j % 3 == 0 && j > 0)
				cell ^= xor_mask;
			ppOutput[i][j] = cell;
		}
	}
	return ppOutput;
}
unsigned long long get_cpu_cycles(morpho_func_t morpho, uint8 **ppInput,  long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
    unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	morpho(ppInput, nrl, nrh, ncl, nch, s, ppOutput);
	end = __rdtsc();
	return end - begin;
}
double **init_benchmark_results(long nb_funcs, long size)
{
	double **results;
	results = (double **) malloc(sizeof(double *) * nb_funcs);
	if (!results)    
		exit_on_error("malloc failed");

	results[0] = (double *) malloc(sizeof(double) * nb_funcs * size);
	if (!results[0]) 
		exit_on_error("malloc failed");
	for (long i = 1; i < nb_funcs; i++) results[i] = results[i - 1] + size;
	return results;

}
void free_benchmark_results(double **results, long nb_funcs)
{
	free((char *)results[0]);
	free((char **)results);
}
uint8 **rand_ui8matrix(long size, p_struct_elem_dim s)
{
	uint8** matrix, rand;
	long nrl, nrh, ncl, nch;
	matrix = ui8matrix(0 + s->nrl, size + s->nrh, 0 + s->ncl, size + s->nch);
	// printf("%ld %ld %ld %ld\n",  0 + s->nrl, size + s->nrh, 0 + s->ncl, size + s->nch);
	for (long i = 0; i < size; i++)
		for (long j = 0; j < size; j++) {
			matrix[i][j] =  ui8rand() % 2;
		}

	return matrix;
}
unsigned long long get_min_cycle(unsigned long long *cycles, long packet_size)
{
	unsigned long long min_cycle;	

	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	min_cycle = cycles[0];

	for (long i = 1; i < packet_size; i++) 
		if (min_cycle > cycles[i]) 
			min_cycle = cycles[i];
	return min_cycle;
}

double **benchmark(morpho_func_t morphos[], p_struct_elem_dim s,long nb_funcs, const long nb_tests, long min_size, long max_size, long step, const char *filename_prefix, int save_output)
{	
	// const long nb_tests    = 100;
	const long packet_size = 3;
	unsigned long long  *cycles;
	unsigned long long   min_cycle_sum;
	uint8 **ppInput, **ppOutput;
	double **results;

	long i = 0, size = 0, k =0, l = 0, cnt = 0;
	long nrl, nrh, ncl, nch;
	long nb_elements = 0;
	char filename[128];


	results = init_benchmark_results(nb_funcs, (max_size - min_size)  + 1);
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);

	
	cnt = 0;
	for (size = min_size; size < max_size + 1; size += step) {
		ppInput = ui8matrix_test_input(size, 1, s);
		ppOutput = ui8matrix(0, size - 1, 0, size - 1);
		ncl = nrl = 0; 
		nch = nrh = size - 1;
		nb_elements = (size * size);
		for (i = 0; i < nb_funcs; i++)  {
			min_cycle_sum = 0;
			for (k = 0; k < nb_tests; k++) {

				for (l = 0; l < packet_size; l++) 
					cycles[l] = get_cpu_cycles(morphos[i], ppInput, nrl, nrh, ncl, nch, s, ppOutput);
				min_cycle_sum += get_min_cycle(cycles, packet_size);
			} 
			results[i][cnt] = (double)(min_cycle_sum / (nb_tests * nb_elements)) ;
			// printf("%ld (%ld) : %1.2lf\n", size, nb_elements, results[i][cnt]);

			if (save_output) {
				printf("%p", filename);
				getchar();
				binary_to_octal_ui8matrix(ppOutput, 0, size - 1, 0, size - 1);
				snprintf(filename, 128, "%s_%ldx%ld.pgm", filename_prefix, size, size);
				SavePGM_ui8matrix(ppOutput, 0, size - 1, 0, size - 1, filename);
			}
		}
		free_ui8matrix(ppOutput, 0         , size - 1     , 0         , size - 1);
		free_ui8matrix(ppInput , 0 + s->nrl, size + s->nrh, 0 + s->ncl, size + s->nch);
		cnt++;
	}
	free(cycles);
	return results;
}

#define PERM_TO_ROW(perm, irow, ncol) ((perm >> irow * ncol) & 0x1FF) 
#define  ROW_TO_COL(row, icol)        ((row >> (icol)) & 0x1)

void test_morpho(morpho_func_t erosion, morpho_func_t dilation)
{
	const long TEST_SE_NROW     = 5L;
	const long TEST_SE_NCOL     = 5L;
	const long TEST_SE_ORIGIN_X = 2L;
	const long TEST_SE_ORIGIN_Y = 2L;
	const long nrl =  TEST_SE_ORIGIN_Y - (TEST_SE_NROW - 1);
	const long nrh = -TEST_SE_ORIGIN_Y + (TEST_SE_NROW - 1);
	const long ncl =  TEST_SE_ORIGIN_X - (TEST_SE_NCOL - 1);
	const long nch = -TEST_SE_ORIGIN_X + (TEST_SE_NCOL - 1);

	p_struct_elem_dim s = compute_struct_elem_dim(TEST_SE_ORIGIN_X, TEST_SE_ORIGIN_Y, 
	 											  TEST_SE_NROW    , TEST_SE_NCOL    );
											   
	puts("Checking the structuring element used for unit tests.");
	assert(s->nrow == 5 /* Structuring Element : Height  = 5 */);
	assert(s->ncol == 5 /* Structuring Element : Width   = 5 */);
	assert(s->y0 == 2 /* Structuring Element : originY = 2 */);
	assert(s->x0 == 2 /* Structuring Element : originX = 2 */);
	assert(nrl == -2    /* Test input rect : nrl == -2*/);
	assert(nrh ==  2    /* Test input rect : nrh ==  2*/);
	assert(ncl == -2    /* Test input rect : ncl == -2*/);
	assert(nch ==  2    /* Test input rect : nch ==  2*/);
	puts("The structuring element has the correct dimension for tests.");
	// uint8 rect[TEST_SE_NROW][TEST_SE_NCOL];

	uint8 **ppInput = ui8matrix(nrl, nrh, ncl, nch);
	
	// A binary 5x5 rectangle has 2^25-1 combinations
	// [0 ~ 1FFFFF]
	uint32 min = 0, max = 0x1FFFFFF; // 2^25
	uint32 perm = 0, row = 0;
	const uint32 print_cnt = 0xFFFFF;

	/**
	 * Convert the permutation to matrix format.
	 * 
	 * For example: 
	 * perm = 0b1111100000111110000011111
	 * =>
	 * ppInput = {1,1,1,1,1,
	 *            0,0,0,0,0,
	 *            1,1,1,1,1,
	 *            0,0,0,0,0,
	 *            1,1,1,1,1};
	 **/

	// Loop 1: Get permutations :
	for (perm = 0; perm < max + 1; perm++) {
		// Loop 2 ~ 3 : Tranform the binary permutation to 5x5 rectangle.
		for (int i = nrl; i < nrh + 1; i++) {
			// Get the i-th row from the permutation.
			row = PERM_TO_ROW(perm, (i - nrl), TEST_SE_NCOL); 
			for (int j = ncl; j < nch + 1; j++) {
				// Get the j-th column from the row.
				ppInput[i][j] = ROW_TO_COL(row, (j - ncl));
			}
		}

		if (perm % print_cnt == 0) 
			printf("%8dth test finished (%08d / %8d).\n", perm, perm, max);
		
		if (perm == max)	assert(Morpho_Test_5x5_Rect( erosion, ppInput, s) == 1 /*perm == 0x1FFFFFF */);
		else				assert(Morpho_Test_5x5_Rect( erosion, ppInput, s) == 0 /*perm <  0x1FFFFFF */);
		if (perm == min)	assert(Morpho_Test_5x5_Rect(dilation, ppInput, s) == 0 /*perm == 0x0000000 */);
		else 				assert(Morpho_Test_5x5_Rect(dilation, ppInput, s) == 1 /*perm >  0x0000000 */);
	}
	puts("Morpho : Passed all tests.");
	free_ui8matrix(ppInput, nrl, nrh, ncl, nch);
	free_structuring_element(s);
}
uint8 Morpho_Test_5x5_Rect(morpho_func_t morpho, uint8 **ppInput, p_struct_elem_dim s)
{
	uint8     output = -1;
	uint8  * pOutput = &output;
	uint8 **ppOutput = &pOutput;
	morpho(ppInput, 0, 0, 0, 0, s, ppOutput);
	return output;
}

// uint8 Erosion_Test_5x5_Rect(uint8 **ppRect, p_struct_elem_dim s)
// {
// 	uint8     output = -1;
// 	uint8  * pOutput = &output;
// 	uint8 **ppOutput = &pOutput;
// 	ui8matrix_erosion_naive(ppRect, 0, 0, 0, 0, s, ppOutput);
// 	return output;
// }

// uint8 Dilation_Test_5x5_Rect(uint8 **ppRect, p_struct_elem_dim s)
// {
// 	uint8     output = -1;
// 	uint8  * pOutput = &output;
// 	uint8 **ppOutput = &pOutput;
// 	ui8matrix_dilation_naive(ppRect, 0, 0, 0, 0, s, ppOutput);
// 	return output;
// }