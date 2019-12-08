#include <nrdef.h>
#include <nrutil.h>
#include <img.h>
#include <morpho.h>
#include <stdbool.h>
#include <test_morpho.h>
#include <stdio.h>
#include <stdlib.h>
#include <util.h>
#include <mutil.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>
#include <time.h>
#include <math.h>
#include <verif_data.h>


unsigned long long get_cpu_cycles(struct morpho_set *ptr_mset, uint8 **ppInput,  long nrl, long nrh, long ncl, long nch, uint8 **ppOutput)
{
    unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	ptr_mset->morpho_func(ppInput, nrl, nrh, ncl, nch, ptr_mset->s, ppOutput);
	end = __rdtsc();
	return end - begin;
}

unsigned long long get_min_cpu_cycles(struct morpho_set *ptr_mset, long packet_size, uint8 **ppInput,  long nrl, long nrh, long ncl, long nch, uint8 **ppOutput)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++) {
		cycles[i] = get_cpu_cycles(ptr_mset, ppInput, nrl, nrh, ncl, nch, ppOutput);
	}

	min_cycles = cycles[0];
	for (long i = 1; i < packet_size; i++) 
		if (min_cycles > cycles[i]) 
			min_cycles = cycles[i];
	free(cycles);
	return min_cycles;
}



double **benchmark(struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests)
{
	const int packet_size = 3;
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	double **results;
	uint8 **ppInput, **ppOutput;
	
	nrl = ncl = 0;
	nrh = nch = hs - 1;
	// Find the largest margins.
	for (idx_set = 0; idx_set < nb_sets; idx_set++) {
		nrl = min(msets[idx_set].s->nrl, nrl);
		ncl = min(msets[idx_set].s->ncl, ncl);
		nrh = max(msets[idx_set].s->nrh + hs - 1, nrh);
		nch = max(msets[idx_set].s->nch + hs - 1, nch);
	}
	
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1);
	// the least largest 3x3 checkered square matrix as an input / output matrix.
	ppInput  = ui8matrix_checker(nrl, nrh, ncl, nch, 3, 1); 
	ppOutput = ui8matrix(nrl, nrh, ncl, nch);
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		for (idx_set = 0; idx_set < nb_sets; idx_set++) {
			begin = __rdtsc();			
			memset_ui8matrix(ppOutput, 0, 0, size, 0, size);

			min_cycles_sum = 0;
			for (idx_test = 0; idx_test < nb_tests; idx_test++)
				min_cycles_sum += get_min_cpu_cycles(&msets[idx_set], packet_size, ppInput, 0, size, 0, size, ppOutput);
			
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * packet_size * (size + 1) * (size + 1)));
			end = __rdtsc();
			printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles.\n",  msets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin) / (nb_tests * 3));			
		}

		cnt++;
	}
	free_ui8matrix(ppInput, nrl, nrh, ncl, nch);
	free_ui8matrix(ppOutput,nrl, nrh, ncl, nch);
	return results;
}


bool check_dimension_of_square_structuring_element(p_struct_elem_dim s,  long size)
{
	const long SE_SIZE = size;
	const long SE_ORIGIN = size / 2;
	const long nrl =  SE_ORIGIN - (SE_SIZE - 1);
	const long nrh = -SE_ORIGIN + (SE_SIZE - 1);
	const long ncl =  SE_ORIGIN - (SE_SIZE - 1);
	const long nch = -SE_ORIGIN + (SE_SIZE - 1);
	if(s->nrow != size 		) return false;
	if(s->ncol != size 		) return false;
	if(s->y0   != SE_ORIGIN ) return false;
	if(s->x0   != SE_ORIGIN ) return false;
	if(s->nrl  != nrl       ) return false;
	if(s->nrh  != nrh       ) return false;
	if(s->ncl  != ncl       ) return false;
	if(s->nch  != nch       ) return false;
	return true;
}


#define STRUCTURING_ELEMENT_DIM(s) s->nrl, s->nrh, s->ncl, s->nch
#define PROGRESS_FACTOR 		   10
void print_progress(uint32 current, uint32 max)
{	
	static const int n = 9;
	printf("\tTest progress : [%*d / %-*d].\n", n, current, n, max);
}

void test_implementation_erosion_3x3(struct morpho_set *erosion_set)
{
	assert(check_for_3x3_structuring_element(erosion_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max;
	uint8 **ppInput;
	
	ppInput = ui8matrix(STRUCTURING_ELEMENT_DIM(erosion_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	// Loop 1: Get permutations :
	printf("Implementation test : %s\n", erosion_set->func_name);
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(ppInput, STRUCTURING_ELEMENT_DIM(erosion_set->s), perm);

		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == max) assert(morpho_produces_one(erosion_set, ppInput) == true);
		else			 assert(morpho_produces_one(erosion_set, ppInput) == false);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(ppInput, STRUCTURING_ELEMENT_DIM(erosion_set->s));
}
void test_implementation_dilation_3x3(struct morpho_set *dilation_set)
{
	assert(check_for_3x3_structuring_element(dilation_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max, min = 0;
	uint8 **ppInput;
	
	ppInput = ui8matrix(STRUCTURING_ELEMENT_DIM(dilation_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	printf("Implementation test : %s\n", dilation_set->func_name);
	// Loop 1: Get permutations :
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(ppInput, STRUCTURING_ELEMENT_DIM(dilation_set->s), perm);
		
		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == min) assert(morpho_produces_one(dilation_set, ppInput) == false);
		else			 assert(morpho_produces_one(dilation_set, ppInput) == true);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(ppInput, STRUCTURING_ELEMENT_DIM(dilation_set->s));
}

void test_implementation_erosion_5x5(struct morpho_set *erosion_set)
{
	assert(check_for_5x5_structuring_element(erosion_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 5;
	uint32 perm, max;
	uint8 **ppInput;
	
	ppInput = ui8matrix(STRUCTURING_ELEMENT_DIM(erosion_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	printf("Implementation test : %s\n", erosion_set->func_name);
	// Loop 1: Get permutations :
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(ppInput, STRUCTURING_ELEMENT_DIM(erosion_set->s), perm);

		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == max) assert(morpho_produces_one(erosion_set, ppInput) == true);
		else			 assert(morpho_produces_one(erosion_set, ppInput) == false);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(ppInput, STRUCTURING_ELEMENT_DIM(erosion_set->s));

}
void test_implementation_dilation_5x5(struct morpho_set *dilation_set)
{
	assert(check_for_5x5_structuring_element(dilation_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 5;
	uint32 perm, max, min = 0;
	uint8 **ppInput;
	
	ppInput = ui8matrix(STRUCTURING_ELEMENT_DIM(dilation_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	// Loop 1: Get permutations :
	printf("Implementation test : %s\n", dilation_set->func_name);
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(ppInput, STRUCTURING_ELEMENT_DIM(dilation_set->s), perm);
		
		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == min) assert(morpho_produces_one(dilation_set, ppInput) == false);
		else			 assert(morpho_produces_one(dilation_set, ppInput) == true);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(ppInput, STRUCTURING_ELEMENT_DIM(dilation_set->s));
}
bool morpho_produces_one(struct morpho_set *dilation_set, uint8** ppInput)
{
	uint8     output = -1;
	uint8  * pOutput = &output;
	uint8 **ppOutput = &pOutput;

	dilation_set->morpho_func(ppInput, 0, 0, 0, 0, dilation_set->s, ppOutput);
	return output == 1;
}

uint8** prologue_test_integration_3x3(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput)
{
	long nrl = mset->s->nrl;
	long nrh = mset->s->nrh + NROW_EVAL_OUT_SET_3X3 - 1;
	long ncl = mset->s->ncl;
	long nch = mset->s->nch + NCOL_EVAL_OUT_SET_3X3 - 1;

	ppInput = ui8matrix(nrl, nrh, 
					    ncl, nch);

	ppOutput = ui8matrix(0, NROW_EVAL_OUT_SET_3X3 - 1, 
						 0, NCOL_EVAL_OUT_SET_3X3 - 1);

	memset_ui8matrix(ppInput, 0, nrl, nrh, 
								 ncl, nch);

	memcpy(&ppInput[nrl][ncl], morpho3x3_test_input, SIZE_EVAL_IN_SET_3X3 * sizeof(ppInput[nrl][ncl]));
	printf("Integration test : %s\n", mset->func_name);
	return ppInput;
}

uint8** prologue_test_integration_5x5(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput)
{
	long nrl, nrh, ncl, nch;
	
	nrl = mset->s->nrl; 
	nrh = mset->s->nrh + NROW_EVAL_OUT_SET_5X5 - 1; 
	ncl = mset->s->ncl; 
	nch = mset->s->nch + NCOL_EVAL_OUT_SET_5X5 - 1;
	
	ppInput = ui8matrix(nrl, nrh, 
						ncl, nch);

	memset_ui8matrix(ppInput, 0, nrl, nrh, 
								 ncl, nch);

	memcpy(&ppInput[nrl][ncl], morpho5x5_test_input, SIZE_EVAL_IN_SET_5X5 * sizeof(ppInput[nrl][ncl]));
	printf("Integration test : %s\n", mset->func_name);
	return ppInput;
}

void epilogue_test_integration_3x3(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput)
{
	puts("Passed the integration test.\n");
	free_ui8matrix(ppInput , mset->s->nrl, mset->s->nrh + NROW_EVAL_OUT_SET_3X3 - 1, 
					    	 mset->s->ncl, mset->s->nch + NCOL_EVAL_OUT_SET_3X3 - 1);
	free_ui8matrix(ppOutput, 0, NROW_EVAL_OUT_SET_3X3 - 1, 
						     0, NCOL_EVAL_OUT_SET_3X3 - 1);
}

void epilogue_test_integration_5x5(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput)
{
	puts("Passed the integration test.\n");
	free_ui8matrix(ppInput , mset->s->nrl, mset->s->nrh + NROW_EVAL_OUT_SET_5X5 - 1, 
					    	 mset->s->ncl, mset->s->nch + NCOL_EVAL_OUT_SET_5X5 - 1);
	free_ui8matrix(ppOutput, 0, NROW_EVAL_OUT_SET_5X5 - 1, 
						     0, NCOL_EVAL_OUT_SET_5X5 - 1);
}


void test_integration_erosion_3x3(struct morpho_set *erosion_set)
{
	check_dimension_of_square_structuring_element(erosion_set->s, 3);
	uint8 **ppInput, **ppOutput;

	ppInput  = prologue_test_integration_3x3(erosion_set, ppInput, ppOutput);
	ppOutput = ui8matrix(0, NROW_EVAL_OUT_SET_3X3 - 1, 
						 0, NCOL_EVAL_OUT_SET_3X3 - 1);


	erosion_set->morpho_func(ppInput, 0, NROW_EVAL_OUT_SET_3X3 - 1, 
									  0, NCOL_EVAL_OUT_SET_3X3 - 1, 
							 		  erosion_set->s,
							 		  ppOutput);
	assert(!memcmp(ppOutput[0], erosion3x3_test_output, SIZE_EVAL_OUT_SET_3X3 * sizeof(**ppOutput)));

	
						 	  	
	epilogue_test_integration_3x3(erosion_set, ppInput, ppOutput);

}
void test_integration_erosion_5x5(struct morpho_set *erosion_set)
{
	check_dimension_of_square_structuring_element(erosion_set->s, 5);
	uint8 **ppInput, **ppOutput;

	ppInput  = prologue_test_integration_5x5(erosion_set, ppInput, ppOutput);
	ppOutput = ui8matrix(0, NROW_EVAL_OUT_SET_5X5 - 1, 
						 0, NCOL_EVAL_OUT_SET_5X5 - 1);

	erosion_set->morpho_func(ppInput, 0, NROW_EVAL_OUT_SET_5X5 - 1, 
									  0, NCOL_EVAL_OUT_SET_5X5 - 1, 
							 		  erosion_set->s,
							 		  ppOutput);
									   
	assert(!memcmp(ppOutput[0], erosion5x5_test_output, SIZE_EVAL_OUT_SET_5X5 * sizeof(**ppOutput)));
	epilogue_test_integration_5x5(erosion_set, ppInput, ppOutput);
}
void test_integration_dilation_3x3(struct morpho_set *dilation_set)
{
	check_dimension_of_square_structuring_element(dilation_set->s, 3);
	uint8 **ppInput, **ppOutput;
	
	ppInput  = prologue_test_integration_3x3(dilation_set, ppInput, ppOutput);
	ppOutput = ui8matrix(0, NROW_EVAL_OUT_SET_3X3 - 1, 
						 0, NCOL_EVAL_OUT_SET_3X3 - 1);

	dilation_set->morpho_func(ppInput, 0, NROW_EVAL_OUT_SET_3X3 - 1, 
									   0, NCOL_EVAL_OUT_SET_3X3 - 1, 
							 		   dilation_set->s,
							 		   ppOutput);
	display_ui8matrix(ppOutput, 0, NROW_EVAL_OUT_SET_3X3 - 1, 
						 		0, NCOL_EVAL_OUT_SET_3X3 - 1, "%u", "Check");
	assert(!memcmp(ppOutput[0], dilation3x3_test_output, SIZE_EVAL_OUT_SET_3X3 * sizeof(**ppOutput)));

	epilogue_test_integration_3x3(dilation_set, ppInput, ppOutput);
}
void test_integration_dilation_5x5(struct morpho_set *dilation_set)
{
	check_dimension_of_square_structuring_element(dilation_set->s, 5);
	uint8 **ppInput, **ppOutput;
	
	ppInput  = prologue_test_integration_5x5(dilation_set, ppInput, ppOutput);
	ppOutput = ui8matrix(0, NROW_EVAL_OUT_SET_5X5 - 1, 
						 0, NCOL_EVAL_OUT_SET_5X5 - 1);

	dilation_set->morpho_func(ppInput, 0, NROW_EVAL_OUT_SET_5X5 - 1, 
									   0, NCOL_EVAL_OUT_SET_5X5 - 1, 
							 		   dilation_set->s,
							 		   ppOutput);
	assert(!memcmp(ppOutput[0], dilation5x5_test_output, SIZE_EVAL_OUT_SET_5X5 * sizeof(**ppOutput)));
	epilogue_test_integration_5x5(dilation_set, ppInput, ppOutput);
}


uint8**  ui8matrix_permutation (uint8** m, long nrl, long nrh, long ncl, long nch, uint32 perm)
{
	/**
	 * Convert the permutation to a matrix format.
	 * 
	 * For example: 
	 * perm = 0b111000111
	 * =>
	 * ppInput = {1,1,1,
	 *            0,0,0,
	 *            1,1,1};
	 **/
	uint32 extracted_bits = 0;
	for (int i = nrl; i < nrh + 1; i++) {
		// Extract (i) th ~ (i + nrow) th bits from permutation.
		extracted_bits = extract_bits_from_permutation(perm, (i - nrl), (nch - ncl + 1)); 
		for (int j = ncl; j < nch + 1; j++) {
			// Get the j-th column from the row.
			m[i][j] = get_column_at(extracted_bits, (j - ncl));
		}
	}

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

uint8 **ui8matrix_checker(long nrl, long nrh, long ncl, long nch, const long chkr_size,  const uint8 xor_mask)
{
	uint8** ppOutput, cell;
	long zigzag, row0, col0;

	ppOutput = ui8matrix(nrl, nrh, ncl, nch);
	memset_ui8matrix(ppOutput, xor_mask, nrl, nrh, ncl, nch);
	
	cell = xor_mask;
	
	row0 = (long)ceil(fabs(nrl) / chkr_size) * chkr_size * (nrl < 0 ? -1 : 1);
	col0 = (long)ceil(fabs(ncl) / chkr_size) * chkr_size * (ncl < 0 ? -1 : 1);
	
	for (long row = row0; row < nrh + 1; row+=chkr_size) {
		// change the initial position to make a zigzag pattern.
		zigzag = ((zigzag + 1) % 2 * chkr_size);
		for (long col = col0 + zigzag; col < nch + 1; col+=chkr_size * 2) {
			// draw chkr_size x chkr_size size checker :
			for (long y = 0; y < chkr_size; y++) 
				for (long x = 0; x < chkr_size; x++) {
					// Check bound :
					if ((row + y) >= nrl && (col + x) >= ncl &&
					    (row + y) <= nrh && (col + x) <= nch)
						ppOutput[row + y][col + x] ^= xor_mask;
			}
		}
	}
	return ppOutput;
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
