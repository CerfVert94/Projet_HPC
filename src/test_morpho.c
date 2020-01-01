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
			
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			end = __rdtsc();
			printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles.\n",  msets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin));			
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
bool morpho_produces_one(struct morpho_set *mset, uint8** ppInput)
{
	uint8 X[5][5]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	uint8 *Y[5] = {X[0] + 2, X[1] + 2, X[2] + 2, X[3] + 2, X[4] + 2};
	uint8 **Z = Y + 2;

	mset->morpho_func(ppInput, 0, 0, 0, 0, mset->s, Z);
	return Z[0][0] == 1;
}

uint8** prologue_test_integration_3x3(struct morpho_set *mset, uint8** ppInput, uint8* test_case , long nrl, long nrh, long ncl, long nch, uint8** ppOutput)
{
	const long size = (nrh - nrl + 1) * (nch - ncl + 1);
	

	ppInput = ui8matrix(nrl, nrh, 
					    ncl, nch);

	memset_ui8matrix(ppInput, 0, nrl, nrh, 
								 ncl, nch);

	memcpy(&ppInput[nrl][ncl], test_case, size * sizeof(ppInput[nrl][ncl]));

	// display_ui8matrix(ppInput , nrl, nrh, 
						 		// ncl, nch, "%u", "Input");
	printf("Integration test : %s\n", mset->func_name);
	return ppInput;
}


void epilogue_test_integration_3x3(struct morpho_set *mset, uint8** ppInput, long in_nrl, long in_nrh, long in_ncl, long in_nch,  long out_nrl, long out_nrh, long out_ncl, long out_nch, uint8** ppOutput)
{
	puts("Passed the integration test.\n");
	free_ui8matrix(ppInput , in_nrl, in_nrh, 
					    	 in_ncl, in_nch);
	free_ui8matrix(ppOutput, out_nrl, out_nrh, 
						     out_ncl, out_nch);
}



void test_integration_morpho_3x3(struct morpho_set *mset, uint8** X, long img_nrl, long img_nrh, long img_ncl,long img_nch, uint8** correct_output)
{
    struct struct_elem_dim *s = mset->s;
	char filename[128];
	static int i = 0;
	long nrl = img_nrl + s->nrl,\
		 nrh = img_nrh + s->nrh,\
		 ncl = img_ncl + s->ncl,\
		 nch = img_nch + s->nch;
    long size = (nrh - nrl + 1) * (nch - ncl + 1);
    uint8** Y = ui8matrix(nrl, nrh, ncl, nch);
	memset_ui8matrix(Y, 0, nrl, nrh, ncl, nch);
	
	mset->morpho_func(X, img_nrl, img_nrh, img_ncl, img_nch, mset->s, Y);
	assert(!memcmp(Y[nrl] + ncl, correct_output[nrl] + ncl, size * sizeof(**Y)));
	// snprintf(filename, 128, "output/car[%d][%s].pgm",i++, mset->func_name);
    // SavePGM_ui8matrix(Y, nrl, nrh, ncl, nch, filename);
    free_ui8matrix(Y, nrl, nrh, ncl, nch);
}

void test_erosions  (struct morpho_set *erosion_sets , const int nb_implementations, bool display)
{
    const long start_x = 280, start_y = 200, end_x=320, end_y=240;
    long x = 0, y = 0, nrl, nrh, ncl, nch;
    char filename[128];
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_erosion_naive", ui8matrix_erosion_naive, .s= erosion_sets[0].s};
    for (int i = 0; i < nb_implementations; i++) {
        if (erosion_sets[i].s->ncol == 3) {
            test_implementation_erosion_3x3(&erosion_sets[i]);
        }
        else {//if (erosion_sets[i].s == s5)
            // test_implementation_erosion_5x5(&erosion_sets[i]);
            // test_integration_erosion_5x5(&erosion_sets[i]);
        }
    }
    test_intergration("../car3/car_3000.pgm", &naive_morpho_set, erosion_sets, nb_implementations, display);

  
    // snprintf(filename, 128, "testcases/car{%ld}{%ld}.pgm",x,y);
}
void test_dilations(struct morpho_set *dilation_sets, const int nb_implementations, bool display)
{
    const long start_x = 280, start_y = 200, end_x=320, end_y=240;
    long x = 0, y = 0, nrl, nrh, ncl, nch;
    char filename[128];
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_dilation_naive", ui8matrix_dilation_naive, .s= dilation_sets[0].s};
    for (int i = 0; i < nb_implementations; i++) {
        if (dilation_sets[i].s->ncol == 3) {
            test_implementation_dilation_3x3(&dilation_sets[i]);
        }
        else {//if (dilation_sets[i].s == s5)
            // test_implementation_dilation_5x5(&dilation_sets[i]);
            // test_integration_dilation_5x5(&dilation_sets[i]);
        }
    }
    test_intergration("../car3/car_3000.pgm", &naive_morpho_set, dilation_sets, nb_implementations, display);

  
    // snprintf(filename, 128, "testcases/car{%ld}{%ld}.pgm",x,y);
}
void test_sequences(struct morpho_set *sequence_sets, const int nb_implementations, bool display)
{
    const long start_x = 280, start_y = 200, end_x=320, end_y=240;
    long x = 0, y = 0, nrl, nrh, ncl, nch;
    char filename[128];
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_sequence_naive", ui8matrix_sequence_naive, .s= sequence_sets[0].s};
    test_intergration("../car3/car_3000.pgm", &naive_morpho_set, sequence_sets, nb_implementations, display);
    // snprintf(filename, 128, "testcases/car{%ld}{%ld}.pgm",x,y);
}
void test_intergration(char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, const int nb_implementations, bool display)
{
      
    long nrl, ncl, nrh, nch, temp_nrh, temp_nch;
    struct struct_elem_dim *s = naive_morpho_set->s;
    uint8 **image, **X, **Y, **Z;

    image = LoadPGM_ui8matrix(filename, &nrl, &nrh, &ncl, &nch);

    // Test Input
    X = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
    // Full zero intialization & copy image 
    memset_ui8matrix(X, 0, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
    copy_ui8matrix_ui8matrix(image, nrl, nrh, ncl, nch, X);

    for(temp_nrh = nrl + 9; temp_nrh < nrh + 1; temp_nrh++){
        for(temp_nch = ncl + 9; temp_nch < nch + 1; temp_nch++){
            // Test Output
            Y = ui8matrix(nrl, temp_nrh, ncl, temp_nch);
            // Valid Output
            Z = ui8matrix(nrl, temp_nrh, ncl, temp_nch);
            memset_ui8matrix(Z, 0, nrl, temp_nrh , ncl, temp_nch); 
            naive_morpho_set->morpho_func(X, nrl, temp_nrh, ncl, temp_nch, s, Z);
            if (display) {
                display_ui8matrix(X, nrl, temp_nrh, ncl, temp_nch, "%03u ", "Input");
                display_ui8matrix(Z, nrl, temp_nrh, ncl, temp_nch, "%03u ", naive_morpho_set->func_name);
            }

            for (int i = 0; i < nb_implementations; i++){
                assert(morpho_sets[i].s->nrow == 3 && morpho_sets[i].s->ncol == 3);
                
                printf("Integration test : "LALIGNED_STR" (%-30s)\n", morpho_sets[i].func_name, filename);
                printf("\tTesting for %ld x %ld\n", temp_nch + 1, temp_nrh + 1);

                memset_ui8matrix(Y, 0, nrl, temp_nrh, ncl, temp_nch); 
                morpho_sets[i].morpho_func(X, nrl, temp_nrh, ncl, temp_nch, morpho_sets[i].s, Y);
            
                if (display) {
                    display_ui8matrix(X, nrl, temp_nrh, ncl, temp_nch, "%03u ", "Input");
                    display_ui8matrix(Y, nrl, temp_nrh, ncl, temp_nch, "%03u ", morpho_sets[i].func_name);
                    display_ui8matrix(Z, nrl, temp_nrh, ncl, temp_nch, "%03u ", naive_morpho_set->func_name);
                }
                assert(!memcmp(Y[nrl] + ncl, Z[nrl] + ncl, (temp_nrh - nrl + 1) * (temp_nch - ncl + 1) * sizeof(**Z)));
	            printf("\tTest passed\n");
            }
            free_ui8matrix(Y, nrl, temp_nrh , ncl, temp_nch); 
            free_ui8matrix(Z, nrl, temp_nrh , ncl, temp_nch); 
        }
    }
    free_ui8matrix(X, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
    
}


// void test_integration_erosion_3x3_r1(struct morpho_set *erosion_set)
// {
// 	long in_nrl = NRL_R1EVAL_IN_SET_3X3, in_nrh = NRH_R1EVAL_IN_SET_3X3, in_ncl = NCL_R1EVAL_IN_SET_3X3, in_nch = NCH_R1EVAL_IN_SET_3X3;
// 	long out_nrl = NRL_R1EVAL_OUT_SET_3X3, out_nrh = NRH_R1EVAL_OUT_SET_3X3, out_ncl = NCL_R1EVAL_OUT_SET_3X3, out_nch = NCH_R1EVAL_OUT_SET_3X3;
// 	long out_size = (out_nrh - out_nrl + 1) * (out_nch - out_ncl + 1);
// 	uint8 *test_case = morpho3x3_r1_test_input;
// 	uint8 *correct_output = erosion3x3_r1_test_output;
// 	struct morpho_set *mset = erosion_set;
// 	struct struct_elem_dim *s = mset->s;
// 	uint8 **ppInput, **ppOutput;

// 	check_dimension_of_square_structuring_element(s, 3);
	
// 	ppInput  = prologue_test_integration_3x3(mset, ppInput, test_case, in_nrl, in_nrh, in_ncl, in_nch, ppOutput);
// 	ppOutput = ui8matrix(out_nrl, out_nrh, 
// 						 out_ncl, out_nch);
// 	memset_ui8matrix(ppOutput, 0, in_nrl, in_nrh, 
// 								  in_ncl, in_nch);

// 	mset->morpho_func(ppInput, out_nrl - s->nrl, out_nrh - s->nrh, 
// 				 			   out_ncl - s->ncl, out_nch - s->nch, 
// 				 	 		   mset->s,
// 				 	 		   ppOutput);

// 	assert(!memcmp(ppOutput[out_nrl] + out_ncl, correct_output, out_size * sizeof(**ppOutput)));

// }
// void test_integration_erosion_3x3_r2(struct morpho_set *erosion_set)
// {
// 	long in_nrl = NRL_R2EVAL_IN_SET_3X3, in_nrh = NRH_R2EVAL_IN_SET_3X3, in_ncl = NCL_R2EVAL_IN_SET_3X3, in_nch = NCH_R2EVAL_IN_SET_3X3;
// 	long out_nrl = NRL_R2EVAL_OUT_SET_3X3, out_nrh = NRH_R2EVAL_OUT_SET_3X3, out_ncl = NCL_R2EVAL_OUT_SET_3X3, out_nch = NCH_R2EVAL_OUT_SET_3X3;
// 	long out_size = (out_nrh - out_nrl + 1) * (out_nch - out_ncl + 1);
// 	uint8 *test_case = morpho3x3_r2_test_input;
// 	uint8 *correct_output = erosion3x3_r2_test_output;
// 	struct morpho_set *mset = erosion_set;
// 	struct struct_elem_dim *s = mset->s;
// 	uint8 **ppInput, **ppOutput;

// 	check_dimension_of_square_structuring_element(s, 3);
	
// 	ppInput  = prologue_test_integration_3x3(mset, ppInput, test_case,  in_nrl, in_nrh, in_ncl, in_nch, ppOutput);
// 	ppOutput = ui8matrix(out_nrl, out_nrh, 
// 						 out_ncl, out_nch);
// 	memset_ui8matrix(ppOutput, 0, in_nrl, in_nrh, 
// 								  in_ncl, in_nch);

// 	mset->morpho_func(ppInput, out_nrl - s->nrl, out_nrh - s->nrh, 
// 				 			   out_ncl - s->ncl, out_nch - s->nch, 
// 				 	 		   mset->s,
// 				 	 		   ppOutput);

// 	assert(!memcmp(ppOutput[out_nrl] + out_ncl, correct_output, out_size * sizeof(**ppOutput)));

// 	epilogue_test_integration_3x3(mset, ppInput,  in_nrl,  in_nrh,  in_ncl,  in_nch, 
// 												 out_nrl, out_nrh, out_ncl, out_nch, ppOutput);

// }
// void test_integration_dilation_3x3_r1(struct morpho_set *dilation_set)
// {
// 	long in_nrl = NRL_R1EVAL_IN_SET_3X3, in_nrh = NRH_R1EVAL_IN_SET_3X3, in_ncl = NCL_R1EVAL_IN_SET_3X3, in_nch = NCH_R1EVAL_IN_SET_3X3;
// 	long out_nrl = NRL_R1EVAL_OUT_SET_3X3, out_nrh = NRH_R1EVAL_OUT_SET_3X3, out_ncl = NCL_R1EVAL_OUT_SET_3X3, out_nch = NCH_R1EVAL_OUT_SET_3X3;
// 	long out_size = (out_nrh - out_nrl + 1) * (out_nch - out_ncl + 1);
// 	uint8 *test_case = morpho3x3_r1_test_input;
// 	uint8 *correct_output = dilation3x3_r1_test_output;
// 	struct morpho_set *mset = dilation_set;
// 	struct struct_elem_dim *s = mset->s;
// 	uint8 **ppInput, **ppOutput;

// 	check_dimension_of_square_structuring_element(s, 3);
	
// 	ppInput  = prologue_test_integration_3x3(mset, ppInput, test_case,  in_nrl, in_nrh, in_ncl, in_nch, ppOutput);
// 	ppOutput = ui8matrix(out_nrl, out_nrh, 
// 						 out_ncl, out_nch);
// 	// printf("malloced\n");
// 	memset_ui8matrix(ppOutput, 0, in_nrl, in_nrh, 
// 								  in_ncl, in_nch);

// 	mset->morpho_func(ppInput, out_nrl - s->nrl, out_nrh - s->nrh, 
// 				 			   out_ncl - s->ncl, out_nch - s->nch, 
// 				 	 		   mset->s,
// 				 	 		   ppOutput);

// 	// display_ui8matrix(ppInput, out_nrl, out_nrh, out_ncl, out_nch, "%u", "DR1:");
// 	// display_ui8matrix(ppOutput, out_nrl, out_nrh, out_ncl, out_nch, "%u", "DR1:");
// 	assert(!memcmp(ppOutput[out_nrl] + out_ncl, correct_output, out_size * sizeof(**ppOutput)));

// 	epilogue_test_integration_3x3(mset, ppInput,  in_nrl,  in_nrh,  in_ncl,  in_nch, 
// 												 out_nrl, out_nrh, out_ncl, out_nch, ppOutput);
// }
// void test_integration_dilation_3x3_r2(struct morpho_set *dilation_set)
// {
// 	long in_nrl = NRL_R2EVAL_IN_SET_3X3, in_nrh = NRH_R2EVAL_IN_SET_3X3, in_ncl = NCL_R2EVAL_IN_SET_3X3, in_nch = NCH_R2EVAL_IN_SET_3X3;
// 	long out_nrl = NRL_R2EVAL_OUT_SET_3X3, out_nrh = NRH_R2EVAL_OUT_SET_3X3, out_ncl = NCL_R2EVAL_OUT_SET_3X3, out_nch = NCH_R2EVAL_OUT_SET_3X3;
// 	long out_size = (out_nrh - out_nrl + 1) * (out_nch - out_ncl + 1);
// 	uint8 *test_case = morpho3x3_r2_test_input;
// 	uint8 *correct_output = dilation3x3_r2_test_output;
// 	struct morpho_set *mset = dilation_set;
// 	struct struct_elem_dim *s = mset->s;
// 	uint8 **ppInput, **ppOutput;

// 	check_dimension_of_square_structuring_element(s, 3);
	
// 	ppInput  = prologue_test_integration_3x3(mset, ppInput, test_case,  in_nrl, in_nrh, in_ncl, in_nch, ppOutput);
// 	ppOutput = ui8matrix(out_nrl, out_nrh, 
// 						 out_ncl, out_nch);
// 	memset_ui8matrix(ppOutput, 0, in_nrl, in_nrh, 
// 								  in_ncl, in_nch);

// 	mset->morpho_func(ppInput, out_nrl - s->nrl, out_nrh - s->nrh, 
// 				 			   out_ncl - s->ncl, out_nch - s->nch, 
// 				 	 		   mset->s,
// 				 	 		   ppOutput);

// 	assert(!memcmp(ppOutput[out_nrl] + out_ncl, correct_output, out_size * sizeof(**ppOutput)));
// }
// void epilogue_test_integration_5x5(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput)
// {
// 	puts("Passed the integration test.\n");
// 	free_ui8matrix(ppInput , NRL_EVAL_IN_SET_5X5, NRH_EVAL_IN_SET_5X5, 
// 							 NCL_EVAL_IN_SET_5X5, NCH_EVAL_IN_SET_5X5);
// 	free_ui8matrix(ppOutput, NRL_EVAL_OUT_SET_5X5, NRH_EVAL_OUT_SET_5X5, 
// 							 NCL_EVAL_OUT_SET_5X5, NCH_EVAL_OUT_SET_5X5);
// }
// uint8** prologue_test_integration_5x5(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput)
// {
// 	long nrl = NRL_EVAL_IN_SET_5X5;
// 	long nrh = NRH_EVAL_IN_SET_5X5;
// 	long ncl = NCL_EVAL_IN_SET_5X5;
// 	long nch = NCH_EVAL_IN_SET_5X5;
	
// 	ppInput = ui8matrix(nrl, nrh, 
// 						ncl, nch);

// 	memset_ui8matrix(ppInput, 0, nrl, nrh, 
// 								 ncl, nch);

// 	memcpy(&ppInput[nrl][ncl], morpho5x5_test_input, SIZE_EVAL_IN_SET_5X5 * sizeof(ppInput[nrl][ncl]));
// 	printf("Integration test : %s\n", mset->func_name);
// 	return ppInput;
// }
// void test_integration_erosion_5x5(struct morpho_set *erosion_set)
// {
// 	check_dimension_of_square_structuring_element(erosion_set->s, 5);
// 	uint8 **ppInput, **ppOutput;
// 	struct struct_elem_dim *s = erosion_set->s;
	
// 	ppInput  = prologue_test_integration_5x5(erosion_set, ppInput, ppOutput);

// 	ppOutput = ui8matrix(NRL_EVAL_OUT_SET_5X5, NRH_EVAL_OUT_SET_5X5, 
// 						 NCL_EVAL_OUT_SET_5X5, NCH_EVAL_OUT_SET_5X5);

// 	memset_ui8matrix(ppOutput, 0, NRL_EVAL_OUT_SET_5X5, NRH_EVAL_OUT_SET_5X5, 
// 								  NCL_EVAL_OUT_SET_5X5, NCH_EVAL_OUT_SET_5X5);

// 	erosion_set->morpho_func(ppInput, 0, NRH_EVAL_OUT_SET_5X5 - s->nrh, 
// 									  0, NCH_EVAL_OUT_SET_5X5 - s->nch, 
// 							 		  erosion_set->s,
// 							 		  ppOutput);
									   
// 	assert(!memcmp(ppOutput[NRL_EVAL_OUT_SET_5X5] + NCL_EVAL_OUT_SET_5X5, 
// 				   erosion5x5_test_output, 
// 				   SIZE_EVAL_OUT_SET_5X5 * sizeof(**ppOutput)));

// 	epilogue_test_integration_5x5(erosion_set, ppInput, ppOutput);
// }

// void test_integration_dilation_5x5(struct morpho_set *dilation_set)
// {
// 	check_dimension_of_square_structuring_element(dilation_set->s, 5);
// 	uint8 **ppInput, **ppOutput;
// 	struct struct_elem_dim *s = dilation_set->s;
	
// 	ppInput  = prologue_test_integration_5x5(dilation_set, ppInput, ppOutput);
// 	ppOutput = ui8matrix(NRL_EVAL_OUT_SET_5X5, NRH_EVAL_OUT_SET_5X5, 
// 						 NCL_EVAL_OUT_SET_5X5, NCH_EVAL_OUT_SET_5X5);

// 	memset_ui8matrix(ppOutput, 0, NRL_EVAL_OUT_SET_5X5, NRH_EVAL_OUT_SET_5X5, 
// 								  NCL_EVAL_OUT_SET_5X5, NCH_EVAL_OUT_SET_5X5);
// 	dilation_set->morpho_func(ppInput, 0, NRH_EVAL_OUT_SET_5X5 - s->nrh, 
// 									   0, NCH_EVAL_OUT_SET_5X5 - s->nch, 
// 							 		   dilation_set->s,
// 							 		   ppOutput);
// 	assert(!memcmp(ppOutput[NRL_EVAL_OUT_SET_5X5] + NCL_EVAL_OUT_SET_5X5, 
// 				   dilation5x5_test_output, 
// 				   SIZE_EVAL_OUT_SET_5X5 * sizeof(**ppOutput)));
// 	epilogue_test_integration_5x5(dilation_set, ppInput, ppOutput);
// }


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
