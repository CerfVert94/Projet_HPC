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
#define STRUCTURING_ELEMENT_DIM(s) SE_NRL, SE_NRH, SE_NCL, SE_NCH
#define PROGRESS_FACTOR 		   10
void print_progress(uint32 current, uint32 max)
{	
	static const int n = 9;
	printf("\tTest progress : [%*d / %-*d].\n", n, current, n, max);
}

void test_implementation_erosion3(struct morpho_set *erosion_set)
{
	// assert(check_for3_structuring_element(erosion_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max;
	uint8 **X;
	
	X = ui8matrix(-1, 1, -1, 1);
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
void test_implementation_dilation3(struct morpho_set *dilation_set)
{
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max, min = 0;
	uint8 **X;
	
	X = ui8matrix(-1, 1, -1, 1);
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

void test_implementation_erosion5(struct morpho_set *erosion_set)
{
	// assert(check_for5_structuring_element(erosion_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 5;
	uint32 perm, max;
	uint8 **X;
	
	X = ui8matrix(-2, 2, -2, 2);
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
void test_implementation_dilation5(struct morpho_set *dilation_set)
{
	// assert(check_for_5x5_structuring_element(dilation_set->s) == true);
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 5;
	uint32 perm, max, min = 0;
	uint8 **X;
	
	X = ui8matrix(-2, 2, -2, 2);
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



void test_erosions  (struct morpho_set *erosion_sets , const int nb_implementations, bool display)
{
    const long start_x = 280, start_y = 200, end_x=320, end_y=240;
    long x = 0, y = 0, nrl, nrh, ncl, nch;
    char filename[128];
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_erosion_naive", ui8matrix_erosion_naive};
    for (int i = 0; i < nb_implementations; i++) {
            test_implementation_erosion3(&erosion_sets[i]);
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
    for (int i = 0; i < nb_implementations; i++) {
        if (dilation_sets[i].pack_type == NO_PACK)
			test_implementation_dilation3(&dilation_sets[i]);
		
        
    }
    test_intergration("../car3/car_3000.pgm", &naive_morpho_set, dilation_sets, nb_implementations, display);
	// test_packed_intergration("../car3/car_3000.pgm", &naive_morpho_set, dilation_sets, nb_implementations, display);
  
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
					memset_ui8matrix          (Y, 0, nrl, temp_nrh, ncl, temp_nch); 
					morpho_sets[i].morpho_func(X   , nrl, temp_nrh, ncl, temp_nch, temp_buffer, Y);
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
