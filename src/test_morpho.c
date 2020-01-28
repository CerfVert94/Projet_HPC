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



#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1

unsigned long long get_cpu_cycles(struct morpho_set *ptr_mset, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
    unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	ptr_mset->morpho_func(X, nrl, nrh, ncl, nch, temp_buffer, Y);
	end = __rdtsc();
	return end - begin;
}

unsigned long long get_min_cpu_cycles(struct morpho_set *ptr_mset, long packet_size, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++) {
		cycles[i] = get_cpu_cycles(ptr_mset, X, nrl, nrh, ncl, nch, temp_buffer, Y);
	}

	min_cycles = cycles[0];
	for (long i = 1; i < packet_size; i++) 
		if (min_cycles > cycles[i]) 
			min_cycles = cycles[i];
	free(cycles);
	return min_cycles;
}



double **benchmark(struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	double **results;
	uint8 **X, **Y, **temp_buffer;
	
	nrl = ncl = 0;
	nrh = nch = hs - 1;
	// Find the largest margins.
	for (idx_set = 0; idx_set < nb_sets; idx_set++) {
		nrl = min(SE_NRL, nrl);
		ncl = min(SE_NCL, ncl);
		nrh = max(SE_NRH + hs - 1, nrh);
		nch = max(SE_NCH + hs - 1, nch);
	}
	
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	// the least largest 3x3 checkered square matrix as an input / output matrix.
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		X            = ui8matrix_checker(-2, size + 2, -2, size + 2, 3, 1); 
		Y            = ui8matrix(0, size, 0, size);
		temp_buffer  = ui8matrix(-2, size + 2, -2, size + 2); 
		
		for (idx_set = 0; idx_set < nb_sets; idx_set++) {
			begin = __rdtsc();			
			memset_ui8matrix(Y, 0, 0, size, 0, size);

			min_cycles_sum = 0;
			for (idx_test = 0; idx_test < nb_tests; idx_test++)
				min_cycles_sum += get_min_cpu_cycles(&msets[idx_set], packet_size, X, 0, size, 0, size, temp_buffer, Y);
			
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			end = __rdtsc();
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles.\n",  msets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin));			
		}
	
		free_ui8matrix(temp_buffer, -2, size + 2, -2, size + 2);
		free_ui8matrix(X, -2, size + 2, -2, size + 2);
		free_ui8matrix(Y, 0, size, 0, size);

		cnt++;
	}
	return results;
}

double **benchmark_compression(struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	long packed_nrl, packed_nrh, packed_ncl, packed_nch, bord;
	double **results;
	uint8 **X, **packedX, **Y, **temp_buffer, **Z;
	
	nrl = ncl = 0;
	nrh = nch = hs - 1;
	// Find the largest margins.
	for (idx_set = 0; idx_set < nb_sets; idx_set++) {
		nrl = min(SE_NRL, nrl);
		ncl = min(SE_NCL, ncl);
		nrh = max(SE_NRH + hs - 1, nrh);
		nch = max(SE_NCH + hs - 1, nch);
	}
	
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	// the least largest 3x3 checkered square matrix as an input / output matrix.
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		packedX      = fcpacked_ui8matrix(0, size, 0, size, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord); 
		temp_buffer  = fcpacked_ui8matrix(0, size, 0, size, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord); 
		Y            = fcpacked_ui8matrix(0, size, 0, size, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord); 
		X            = ui8matrix_checker(-bord, size + bord, -bord, size + bord, 3, 1); 
		Z            = ui8matrix(0, size, 0, size); 
		
		
		for (idx_set = 0; idx_set < nb_sets; idx_set++) {
			memset_ui8matrix(Y, 0, packed_nrl, packed_nrh, packed_ncl, packed_nch);

			begin = __rdtsc();			
			fcpack_ui8matrix_ui8matrix(X, 0, size, 0, size, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord, packedX);
			min_cycles_sum = 0;
			for (idx_test = 0; idx_test < nb_tests; idx_test++)
				min_cycles_sum += get_min_cpu_cycles(&msets[idx_set], packet_size, packedX, packed_nrl, packed_nrh, packed_ncl, packed_nch, temp_buffer, Y);
			
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			unfcpack_ui8matrix_ui8matrix(Y, 0, size, 0, size, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord, Z);
			end = __rdtsc();
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles.\n",  msets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin));			
		}
	
		// free_ui8matrix(temp_buffer, -2, size + 2, -2, size + 2);
		free_packed_ui8matrix(packedX    , packed_nrl, packed_nrh, packed_ncl, packed_nch, bord);
		free_packed_ui8matrix(temp_buffer, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord);
		free_packed_ui8matrix(Y          , packed_nrl, packed_nrh, packed_ncl, packed_nch, bord);
		free_ui8matrix(X, -2, size + 2, -2, size + 2);
		free_ui8matrix(Z, 0, size, 0, size);

		cnt++;
	}
	return results;
}

/*
bool check_dimension_of_square_structuring_element(p_struct_elem_dim  long size)
{
	const long SE_SIZE = size;
	const long SE_ORIGIN = size / 2;
	const long nrl =  SE_ORIGIN - (SE_SIZE - 1);
	const long nrh = -SE_ORIGIN + (SE_SIZE - 1);
	const long ncl =  SE_ORIGIN - (SE_SIZE - 1);
	const long nch = -SE_ORIGIN + (SE_SIZE - 1);
	// if(s->nrow != size 		) return false;
	// if(s->ncol != size 		) return false;
	// if(s->y0   != SE_ORIGIN ) return false;
	// if(s->x0   != SE_ORIGIN ) return false;
	if(SE_NRL  != nrl       ) return false;
	if(SE_NRH  != nrh       ) return false;
	if(SE_NCL  != ncl       ) return false;
	if(SE_NCH  != nch       ) return false;
	return true;
}*/


#define STRUCTURING_ELEMENT_DIM(s) SE_NRL, SE_NRH, SE_NCL, SE_NCH
#define PROGRESS_FACTOR 		   10
void print_progress(uint32 current, uint32 max)
{	
	static const int n = 9;
	printf("\tTest progress : [%*d / %-*d].\n", n, current, n, max);
}

void test_implementation_erosion_3x3(struct morpho_set *erosion_set)
{
	// assert(check_for_3x3_structuring_element(erosion_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max;
	uint8 **X;
	
	X = ui8matrix(STRUCTURING_ELEMENT_DIM(erosion_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	// Loop 1: Get permutations :
	printf("Implementation test : %s\n", erosion_set->func_name);
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(X, STRUCTURING_ELEMENT_DIM(erosion_set->s), perm);

		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == max) assert(morpho_produces_one(erosion_set, X) == true);
		else			 assert(morpho_produces_one(erosion_set, X) == false);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(X, STRUCTURING_ELEMENT_DIM(erosion_set->s));
}
void test_implementation_dilation_3x3(struct morpho_set *dilation_set)
{
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max, min = 0;
	uint8 **X;
	
	X = ui8matrix(STRUCTURING_ELEMENT_DIM(dilation_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	printf("Implementation test : %s\n", dilation_set->func_name);
	// Loop 1: Get permutations :
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(X, STRUCTURING_ELEMENT_DIM(dilation_set->s), perm);
		
		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == min) assert(morpho_produces_one(dilation_set, X) == false);
		else			 assert(morpho_produces_one(dilation_set, X) == true);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(X, STRUCTURING_ELEMENT_DIM(dilation_set->s));
}

void test_implementation_erosion_5x5(struct morpho_set *erosion_set)
{
	// assert(check_for_5x5_structuring_element(erosion_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 5;
	uint32 perm, max;
	uint8 **X;
	
	X = ui8matrix(STRUCTURING_ELEMENT_DIM(erosion_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	printf("Implementation test : %s\n", erosion_set->func_name);
	// Loop 1: Get permutations :
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(X, STRUCTURING_ELEMENT_DIM(erosion_set->s), perm);

		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == max) assert(morpho_produces_one(erosion_set, X) == true);
		else			 assert(morpho_produces_one(erosion_set, X) == false);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(X, STRUCTURING_ELEMENT_DIM(erosion_set->s));

}
void test_implementation_dilation_5x5(struct morpho_set *dilation_set)
{
	// assert(check_for_5x5_structuring_element(dilation_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 5;
	uint32 perm, max, min = 0;
	uint8 **X;
	
	X = ui8matrix(STRUCTURING_ELEMENT_DIM(dilation_set->s));
	max = (1 << (size * size)) - 1; // 2^25

	// Loop 1: Get permutations :
	printf("Implementation test : %s\n", dilation_set->func_name);
	for (perm = 0; perm < max + 1; perm++) {
		ui8matrix_permutation(X, STRUCTURING_ELEMENT_DIM(dilation_set->s), perm);
		
		if (perm % (max / PROGRESS_FACTOR) == 0)
			print_progress(perm, max);

		if (perm == min) assert(morpho_produces_one(dilation_set, X) == false);
		else			 assert(morpho_produces_one(dilation_set, X) == true);
	}
	puts("Passed all tests.\n");
	free_ui8matrix(X, STRUCTURING_ELEMENT_DIM(dilation_set->s));
}
bool morpho_produces_one(struct morpho_set *mset, uint8** W)
{
	uint8 X[5][5]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	uint8 tempX[5][5]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};

	uint8 *Y[5] = {X[0] + 2, X[1] + 2, X[2] + 2, X[3] + 2, X[4] + 2};
	uint8 *tempY[5] = {tempX[0] + 2, tempX[1] + 2, tempX[2] + 2, tempX[3] + 2, tempX[4] + 2};
	uint8 **Z = Y + 2;
	uint8 **tempZ = tempY + 2;

	mset->morpho_func(W, 0, 0, 0, 0, tempZ, Z);
	return Z[0][0] == 1;
}

uint8** prologue_test_integration_3x3(struct morpho_set *mset, uint8** X, uint8* test_case , long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8** Y)
{
	const long size = (nrh - nrl + 1) * (nch - ncl + 1);
	

	X = ui8matrix(nrl, nrh, 
					    ncl, nch);

	memset_ui8matrix(X, 0, nrl, nrh, 
								 ncl, nch);

	memcpy(&X[nrl][ncl], test_case, size * sizeof(X[nrl][ncl]));

	// display_ui8matrix(X , nrl, nrh, 
						 		// ncl, nch, "%u", "Input");
	printf("Integration test : %s\n", mset->func_name);
	return X;
}


void epilogue_test_integration_3x3(struct morpho_set *mset, uint8** X, long in_nrl, long in_nrh, long in_ncl, long in_nch,  long out_nrl, long out_nrh, long out_ncl, long out_nch, uint8** temp_buffer, uint8** Y)
{
	puts("Passed the integration test.\n");
	free_ui8matrix(X , in_nrl,  in_nrh , in_ncl,  in_nch);
	free_ui8matrix(Y, out_nrl, out_nrh, out_ncl, out_nch);
}



void test_integration_morpho_3x3(struct morpho_set *mset, uint8** X, long img_nrl, long img_nrh, long img_ncl,long img_nch, uint8 **temp_buffer, uint8** correct_output)
{
	char filename[128];
	static int i = 0;
	long nrl = img_nrl + SE_NRL,\
		 nrh = img_nrh + SE_NRH,\
		 ncl = img_ncl + SE_NCL,\
		 nch = img_nch + SE_NCH;
    long size = (nrh - nrl + 1) * (nch - ncl + 1);
    uint8** Y = ui8matrix(nrl, nrh, ncl, nch);
	memset_ui8matrix(Y, 0, nrl, nrh, ncl, nch);
	
	mset->morpho_func(X, img_nrl, img_nrh, img_ncl, img_nch, temp_buffer,Y);
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
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_erosion_naive", ui8matrix_erosion_naive};
    for (int i = 0; i < nb_implementations; i++) {
            test_implementation_erosion_3x3(&erosion_sets[i]);
    }
    test_intergration("../car3/car_3000.pgm", &naive_morpho_set, erosion_sets, nb_implementations, display);
	
  
    // snprintf(filename, 128, "testcases/car{%ld}{%ld}.pgm",x,y);
}
void test_dilations(struct morpho_set *dilation_sets, const int nb_implementations, bool display)
{
    const long start_x = 280, start_y = 200, end_x=320, end_y=240;
    long x = 0, y = 0, nrl, nrh, ncl, nch;
    char filename[128];
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_dilation_naive", ui8matrix_dilation_naive};
    // for (int i = 0; i < nb_implementations; i++) {
    //     if (dilation_sets[i].pack_type == NO_PACK)
	// 		test_implementation_dilation_3x3(&dilation_sets[i]);
	// 	// else
	// 	// {
	// 	// 	/* code */
	// 	// }
		
        
    // }
    // test_intergration("../car3/car_3000.pgm", &naive_morpho_set, dilation_sets, nb_implementations, display);
	test_packed_intergration("../car3/car_3000.pgm", &naive_morpho_set, dilation_sets, nb_implementations, display);
  
    // snprintf(filename, 128, "testcases/car{%ld}{%ld}.pgm",x,y);
}
void test_sequences(struct morpho_set *sequence_sets, const int nb_implementations, bool display)
{
    const long start_x = 280, start_y = 200, end_x=0, end_y=240;
    long x = 0, y = 0, nrl, nrh, ncl, nch;
    char filename[128];
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_sequence_naive", ui8matrix_sequence_naive};
    test_intergration("../car3/car_3000.pgm", &naive_morpho_set, sequence_sets, nb_implementations, display);
    // snprintf(filename, 128, "testcases/car{%ld}{%ld}.pgm",x,y);
}

void test_intergration(char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, const int nb_implementations, bool display)
{
      
    long nrl, ncl, nrh, nch, temp_nrh, temp_nch;
	long packed_nrl, packed_ncl, packed_nrh, packed_nch;
    // struct struct_elem_dim *s = naive_morpho_set->s;
    uint8 **image, **X, **Y, **Z,**temp_buffer;
	uint8 **packedX, **packedY, **temp_packed_buffer, **unpacked;

    image = LoadPGM_ui8matrix(filename, &nrl, &nrh, &ncl, &nch);

    // Test Input
    X = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// Middle Buffer
    temp_buffer = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);

    // Full zero intialization & copy image 
    memset_ui8matrix(X, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	memset_ui8matrix(temp_buffer, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	
    copy_ui8matrix_ui8matrix(image, nrl, nrh, ncl, nch, X);
	octal_to_binary_ui8matrix(X, nrl, nrh, ncl, nch);

    for(temp_nrh = nrl + 9; temp_nrh < nrh + 1; temp_nrh++){
        for(temp_nch = ncl + 9; temp_nch < nch + 1; temp_nch++){
            // Test Output
            Y = ui8matrix(nrl, temp_nrh, ncl, temp_nch);
            // Valid Output
            Z = ui8matrix(nrl, temp_nrh, ncl, temp_nch);
            memset_ui8matrix(Z, 0, nrl, temp_nrh , ncl, temp_nch); 
            naive_morpho_set->morpho_func(X, nrl, temp_nrh, ncl, temp_nch, temp_buffer, Z);
         
            for (int i = 0; i < nb_implementations; i++){
                printf("Integration test : "LALIGNED_STR" (%-30s)\n", morpho_sets[i].func_name, filename);
                printf("\tTesting for %ld x %ld\n", temp_nch + 1, temp_nrh + 1);

				if (morpho_sets[i].pack_type == NO_PACK) {
					memset_ui8matrix(Y, 0, nrl, temp_nrh, ncl, temp_nch); 
					morpho_sets[i].morpho_func(X, nrl, temp_nrh, ncl, temp_nch, temp_buffer, Y);
				}
                if (display) {
					printf("%ld %ld %ld %ld\n", nrl, temp_nrh, ncl, temp_nch);
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
	free_ui8matrix(temp_buffer, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
    free_ui8matrix(X, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	free_ui8matrix(image, nrl, nrh, ncl, nch);
}

void test_packed_intergration(char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, const int nb_implementations, bool display)
{
      
    long nrl, ncl, nrh, nch, temp_nrh, temp_nch;
	long packed_nrl, packed_ncl, packed_nrh, packed_nch, bord;
    // struct struct_elem_dim *s = naive_morpho_set->s;
    uint8 **image, **X, **packedX, **packedY, **temp_packed_buffer, **Y, **Z;

    // Test Input
    image = LoadPGM_ui8matrix(filename, &nrl, &nrh, &ncl, &nch);
	octal_to_binary_ui8matrix(image, nrl, nrh, ncl, nch);
    X = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	
    
	
    for(temp_nrh = nrl + 9; temp_nrh < nrh + 1; temp_nrh++){
        for(temp_nch = ncl + 9; temp_nch < nch + 1; temp_nch++){
			memset_ui8matrix(X, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
			copy_ui8matrix_ui8matrix(image, nrl, temp_nrh, ncl, temp_nch, X);

			// Packed Input Matrix
			packedX = fcpacked_ui8matrix(nrl, temp_nrh, ncl, temp_nch, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord);
			// Packed Output Matrix (test output)
			packedY = fcpacked_ui8matrix(nrl, temp_nrh, ncl, temp_nch, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord);
			// Middle Buffer
    		temp_packed_buffer = fcpacked_ui8matrix(nrl, temp_nrh, ncl, temp_nch, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord);
			// Valid Output Matrix
			Z = ui8matrix   (      nrl, temp_nrh, ncl, temp_nch);
			memset_ui8matrix(Z, 0, nrl, temp_nrh, ncl, temp_nch);
			// Unpacked Output Matrix
			Y = ui8matrix   (      nrl, temp_nrh, ncl, temp_nch);
			memset_ui8matrix(Y, 0, nrl, temp_nrh, ncl, temp_nch);
			// Get the valid output
			naive_morpho_set->morpho_func(X, nrl, temp_nrh, ncl, temp_nch, NULL, Z);

			// Pack Input 
			fcpack_ui8matrix_ui8matrix  (image, nrl, temp_nrh, ncl, temp_nch,  packed_nrl,  packed_nrh,  packed_ncl,  packed_nch,  bord, packedX);	

			for (int i = 0; i < nb_implementations; i++){
				printf("Integration test : "LALIGNED_STR" (%-30s)\n", morpho_sets[i].func_name, filename);
                printf("\tTesting for %ld x %ld\n", temp_nch + 1, temp_nrh + 1);

            	memset_ui8matrix(packedY, 0, packed_nrl - bord,  packed_nrh + bord,  packed_ncl - bord,  packed_nch + bord);
				// Get packed output
				morpho_sets[i].morpho_func(packedX,  packed_nrl,  packed_nrh,  packed_ncl,  packed_nch, temp_packed_buffer, packedY);	
				unfcpack_ui8matrix_ui8matrix(packedY, nrl, temp_nrh, ncl, temp_nch, packed_nrl,  packed_nrh,  packed_ncl,  packed_nch, bord, Y);
                if (display) {
					printf("%ld %ld %ld %ld\n", nrl, temp_nrh, ncl, temp_nch);
                    display_ui8matrix(X, nrl, temp_nrh, ncl, temp_nch, "%u ", "Input");
                    display_ui8matrix(Y, nrl, temp_nrh, ncl, temp_nch, "%u ", morpho_sets[i].func_name);
                    display_ui8matrix(Z, nrl, temp_nrh, ncl, temp_nch, "%u ", naive_morpho_set->func_name);
                }
                assert(!memcmp(Y[nrl] + ncl, Z[nrl] + ncl, (temp_nrh - nrl + 1) * (temp_nch - ncl + 1) * sizeof(**Z)));
				printf("\tTest passed\n");
				
			}
            free_packed_ui8matrix(packedX			,  packed_nrl,  packed_nrh,  packed_ncl,  packed_nch,  bord);
			free_packed_ui8matrix(packedY			,  packed_nrl,  packed_nrh,  packed_ncl,  packed_nch,  bord);
			free_packed_ui8matrix(temp_packed_buffer,  packed_nrl,  packed_nrh,  packed_ncl,  packed_nch,  bord);
            free_ui8matrix(Y, nrl, temp_nrh , ncl, temp_nch); 
            free_ui8matrix(Z, nrl, temp_nrh , ncl, temp_nch); 
        }
    }
    free_ui8matrix(X, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	free_ui8matrix(image, nrl, nrh, ncl, nch);
}


uint8**  ui8matrix_permutation (uint8** m, long nrl, long nrh, long ncl, long nch, uint32 perm)
{
	/**
	 * Convert the permutation to a matrix format.
	 * 
	 * For example: 
	 * perm = 0b111000111
	 * =>
	 * X = {1,1,1,
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

double **init_benchmark_results(long nb_funcs, long size, long step)
{
	double **results;
	results = (double **) malloc(sizeof(double *) * nb_funcs);
	if (!results)    
		exit_on_error("malloc failed");

	results[0] = (double *) malloc(sizeof(double) * nb_funcs * (size / step + 1));
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
	uint8** Y, cell;
	long zigzag, row0, col0;

	Y = ui8matrix(nrl, nrh, ncl, nch);
	memset_ui8matrix(Y, xor_mask, nrl, nrh, ncl, nch);
	
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
						Y[row + y][col + x] ^= xor_mask;
			}
		}
	}
	return Y;
}
uint8 **rand_ui8matrix(long size, p_struct_elem_dim s)
{
	uint8** matrix, rand;
	long nrl, nrh, ncl, nch;
	matrix = ui8matrix(0 + SE_NRL, size + SE_NRH, 0 + SE_NCL, size + SE_NCH);
	// printf("%ld %ld %ld %ld\n",  0 + SE_NRL, size + SE_NRH, 0 + SE_NCL, size + SE_NCH);
	for (long i = 0; i < size; i++)
		for (long j = 0; j < size; j++) {
			matrix[i][j] =  ui8rand() % 2;
		}

	return matrix;
}
