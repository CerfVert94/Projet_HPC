#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <nrdef.h>
#include <nrutil.h>
#include <util.h>
#include <img.h>
#include <vnrdef.h>
#include <vnrutil.h>
#include <mynrutil.h>
#include <myvnrutil.h>
#include <img_SIMD.h>
#include <morpho.h>
#include <test_morpho.h>
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

uint8 ui8matrix_erosion3_1pt(uint8**X, long i, long j){
	return X[i - 1][j - 1] & X[i - 1][j + 0] & X[i - 1][j + 1] &
		   X[i + 0][j - 1] & X[i + 0][j + 0] & X[i + 0][j + 1] &
		   X[i + 1][j - 1] & X[i + 1][j + 0] & X[i + 1][j + 1];
}
uint8 ui8matrix_dilation3_1pt(uint8**X, long i, long j){
	return X[i - 1][j - 1] | X[i - 1][j + 0] | X[i - 1][j + 1] |
		   X[i + 0][j - 1] | X[i + 0][j + 0] | X[i + 0][j + 1] |
		   X[i + 1][j - 1] | X[i + 1][j + 0] | X[i + 1][j + 1];
}
uint8 ui8matrix_erosion5_1pt(uint8**X, long i, long j){
	return X[i - 2][j - 2] & X[i - 2][j - 1] & X[i - 2][j + 0] & X[i - 2][j + 1] & X[i - 2][j + 1] & X[i - 2][j + 2] &
		   X[i - 1][j - 2] & X[i - 1][j - 1] & X[i - 1][j + 0] & X[i - 1][j + 1] & X[i - 1][j + 1] & X[i - 1][j + 2] &
		   X[i + 0][j - 2] & X[i + 0][j - 1] & X[i + 0][j + 0] & X[i + 0][j + 1] & X[i + 0][j + 1] & X[i + 0][j + 2] &
		   X[i + 1][j - 2] & X[i + 1][j - 1] & X[i + 1][j + 0] & X[i + 1][j + 1] & X[i + 1][j + 1] & X[i + 1][j + 2] &
		   X[i + 2][j - 2] & X[i + 2][j - 1] & X[i + 2][j + 0] & X[i + 2][j + 1] & X[i + 2][j + 1] & X[i + 2][j + 2];
}
uint8 ui8matrix_dilation5_1pt(uint8**X, long i, long j){
	return X[i - 2][j - 2] | X[i - 2][j - 1] | X[i - 2][j + 0] | X[i - 2][j + 1] | X[i - 2][j + 1] | X[i - 2][j + 2] |
		   X[i - 1][j - 2] | X[i - 1][j - 1] | X[i - 1][j + 0] | X[i - 1][j + 1] | X[i - 1][j + 1] | X[i - 1][j + 2] |
		   X[i + 0][j - 2] | X[i + 0][j - 1] | X[i + 0][j + 0] | X[i + 0][j + 1] | X[i + 0][j + 1] | X[i + 0][j + 2] |
		   X[i + 1][j - 2] | X[i + 1][j - 1] | X[i + 1][j + 0] | X[i + 1][j + 1] | X[i + 1][j + 1] | X[i + 1][j + 2] |
		   X[i + 2][j - 2] | X[i + 2][j - 1] | X[i + 2][j + 0] | X[i + 2][j + 1] | X[i + 2][j + 1] | X[i + 2][j + 2];
}
void test_vec_intergration(uint8 **image, long nrl, long nrh, long ncl, long nch, const char *filename, struct morpho_set *naive_morpho, struct morpho_set *morpho_sets,int nb_sets, bool logging);
void test_vec_intergration(uint8 **image, long nrl, long nrh, long ncl, long nch, const char *filename, struct morpho_set *naive_morpho, struct morpho_set *morpho_sets,int nb_sets, bool logging)
{
      
    long temp_nrh, temp_nch;
    // struct struct_elem_dim *s = naive_morpho_set->s;
	vuint8 **vX, **vY, **vZ,**vTempBuffer;
	vuint8 **vimage;
	vuint8 _mm_equal_ = _mm_set_epi8(0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF);
	uint8 **X, **tempBuffer, **Z;
	int card = card_vuint8();
	int i0, i1, j0, j1;
	assert((nrh - nrl + 1) > 20 && (nch - ncl + 1) > 1);
	// octal_to_binary_ui8matrix(image, nrl, nrh, ncl, nch);
	// scalar image to vector image
	vimage = ui8matrix_to_vui8matrix(image, nrl, nrh, ncl, nch, &i0, &i1, &j0, &j1);
	
	// Test Input
	X		   =  ui8matrix(nrl - 2, nrh + 2, ncl - 2, nch + 2);
	tempBuffer =  ui8matrix(nrl - 2, nrh + 2, ncl - 2, nch + 2);
	Z		   =  ui8matrix(nrl    , nrh    , ncl    , nch    );
    vX 		   = vui8matrix( i0 - 2,  i1 + 2,  j0 - 1,  j1 + 1);
    vY 		   = vui8matrix( i0    ,  i1    ,  j0    ,  j1    );
    // Full zero intialization & copy image 
	zero_vui8matrix(vX, i0 - 2, i1 + 2, j0 - 1, j1 + 1);
	zero_vui8matrix(vY, i0    , i1    , j0    , j1    );
	memset_ui8matrix(X, 0, nrl - 2, nrh + 2, ncl - 2, nch + 2);
	memset_ui8matrix(Z, 0, nrl    , nrh    , ncl    , nch    );
	// zero_vui8matrix(vTempBuffer, nrl - 2, nrh + 2, j0, j1 + 1);
	copy_ui8matrix_ui8matrix(image, nrl, nrh, ncl, nch, X);
    copy_vui8matrix_vui8matrix(vimage, i0, i1, j0, j1, vX);

	naive_morpho->morpho_func(X, nrl, 15, ncl, nch, tempBuffer, Z);
	int v0, v1;
	for (long i = nrh - 20; i < nrh + 1; i++){
		for (long j =  nch - 20; j < nch + 1; j++){
	
	
			naive_morpho->morpho_func(X, nrl, i, ncl, j, tempBuffer, Z);
			s2v1D((int)ncl, (int)j, card, &v0, &v1);
			for (int k = 0; k < nb_sets; k++) {
				printf("Integration test : "LALIGNED_STR" (%-30s)\n", morpho_sets[k].func_name, filename);
				printf("\tTesting for %ld x %ld\n", i + 1, j + 1);

				int l = 0;
				uint8 **Y;
				long nrl_, nrh_, ncl_, nch_;
				morpho_sets->vec_morpho_func(vX, nrl, i, v0, v1 , vY);
				Y = vui8matrix_to_ui8matrix(vY, nrl, i, v0, v1, &nrl_, &nrh_, &ncl_, &nch_);
				
				for (int l = nrl_; l < nrh_ + 1; l++)
					assert(!memcmp_ui8matrix(Z, Y, l, l, ncl_, j));

				printf("Test passed\n");
				free_ui8matrix(Y,nrl_, nrh_, ncl_, nch_);
			}			

		}
	}	
	free_ui8matrix(         X, ncl - 2, nch + 2, nrl - 2, nrh + 2);
	free_ui8matrix(tempBuffer, ncl - 2, nch + 2, nrl - 2, nrh + 2);
	free_ui8matrix(         Z, ncl    , nch    , nrl    , nrh    );
    free_vui8matrix(vX, i0 - 2, i1 + 2, j0 - 1, j1 + 1);
    free_vui8matrix(vY, i0    ,  i1    ,  j0    ,  j1    );
	free_vui8matrix(vimage, i0, i1, j0, j1);
}


void test_erosions(char *filename, struct morpho_set *erosion_sets , const int nb_sets, bool logging)
{
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_erosion_naive", ui8matrix_erosion_naive};
	uint8 **image;
	vuint **vimage;
	long nrl, nrh, ncl, nch;
	int i0, i1, j0, j1;

    image  = LoadPGM_ui8matrix(filename, &nrl, &nrh, &ncl, &nch);
	// vimage  = LoadPGM_vui8matrix(filename, &i0, &i1, &j0, &j1);
		
	for (int i = 0; i < nb_sets; i++) {
        test_implementation_erosion3(&erosion_sets[i]);
	}
	test_intergration(    image, nrl, nrh, ncl, nch, filename, &naive_morpho_set, erosion_sets, nb_sets, logging);
	test_vec_intergration(image, nrl, nrh, ncl, nch, filename, &naive_morpho_set, erosion_sets, nb_sets, logging);
	free_ui8matrix(image, nrl, nrh, ncl, nch);


}
void test_dilations(char *filename, struct morpho_set *dilation_sets, const int nb_sets, bool logging)
{
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_dilation_naive", ui8matrix_dilation_naive};
	uint8 **image;
	long nrl, nrh, ncl, nch;
    image = LoadPGM_ui8matrix(filename, &nrl, &nrh, &ncl, &nch);
    for (int i = 0; i < nb_sets; i++) 
		test_implementation_dilation3(&dilation_sets[i]);
	test_intergration(    image, nrl, nrh, ncl, nch, filename, &naive_morpho_set, dilation_sets, nb_sets, logging);
	test_vec_intergration(image, nrl, nrh, ncl, nch, filename, &naive_morpho_set, dilation_sets, nb_sets, logging);
		
	free_ui8matrix(image, nrl, nrh, ncl, nch);
}
void test_sequences(char *filename, struct morpho_set *sequence_sets, const int nb_sets, bool logging)
{
    struct morpho_set naive_morpho_set = {.func_name = "ui8matrix_sequence_naive", ui8matrix_sequence_naive};
    // test_intergration(filename, &naive_morpho_set, sequence_sets, logging);
}

void print_progress(uint32 current, uint32 max)
{	
	static const int n = 9;
	printf("\tTest progress : [%*d / %-*d].\n", n, current, n, max);
}

void test_implementation_erosion3(struct morpho_set *erosion_set)
{
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max, min = 0;
	uint8 **X;
	vuint8 **vX;
	
	X = ui8matrix(-1, 1, -1, 1);
	vX = vui8matrix(-1, 1, -1, 1);
	max = (1 << (size * size)) - 1; // 2^25

	printf("Implementation test : %s\n", erosion_set->func_name);
	if (erosion_set->instr_type == SCALAR) {
		for (perm = 0; perm < max + 1; perm++) {
			ui8matrix_permutation(X, -1, 1, -1, 1, perm);
			
			if (perm % (max / PROGRESS_FACTOR) == 0)
				print_progress(perm, max);

			if (perm == max) assert(morpho_produces_one(erosion_set, X) == true);
			else			 assert(morpho_produces_one(erosion_set, X) == false);
		}
	}
	else if (erosion_set->instr_type == SIMD) {
		for (perm = 0; perm < max + 1; perm++) {
			vui8matrix_permutation(vX, -1, 1, -1, 1, perm);
			
			if (perm % (max / PROGRESS_FACTOR) == 0)
				print_progress(perm, max);

			
			if (perm == max) {assert(vec_morpho_produces_one(erosion_set, vX) == true); /*printf("True\n");*/}
			else			 {assert(vec_morpho_produces_one(erosion_set, vX) == false);/*printf("False\n");*/}
		}
	}
	puts("Passed all tests.\n");
	free_ui8matrix(X, -1, 1, -1, 1);
	free_vui8matrix(vX, -1, 1, -1, 1);
}
void test_implementation_dilation3(struct morpho_set *dilation_set)
{
	// A binary square matrix has 2^(size*size)-1 combinations
	const int size = 3;
	uint32 perm, max, min = 0;
	uint8 **X;
	vuint8 **vX;
	
	X = ui8matrix(-1, 1, -1, 1);
	vX = vui8matrix(-1, 1, -1, 1);
	max = (1 << (size * size)) - 1; // 2^25

	printf("Implementation test : %s\n", dilation_set->func_name);
	if (dilation_set->instr_type == SCALAR) {
		for (perm = 0; perm < max + 1; perm++) {
			ui8matrix_permutation(X, -1, 1, -1, 1, perm);
			
			if (perm % (max / PROGRESS_FACTOR) == 0)
				print_progress(perm, max);

			if (perm == min) {assert(morpho_produces_one(dilation_set, X) == false); /*printf("False\n");*/}
			else			 {assert(morpho_produces_one(dilation_set, X) == true); /*printf("True\n");*/}
		}
	}
	else if (dilation_set->instr_type == SIMD) {
		for (perm = 0; perm < max + 1; perm++) {
			vui8matrix_permutation(vX, -1, 1, -1, 1, perm);
			
			if (perm % (max / PROGRESS_FACTOR) == 0)
				print_progress(perm, max);

			if (perm == min) {assert(vec_morpho_produces_one(dilation_set, vX) == false);/*printf("False\n");*/}
			else			 {assert(vec_morpho_produces_one(dilation_set, vX) == true);/*printf("True\n");*/}
		}
	}
	puts("Passed all tests.\n");
	free_ui8matrix(X, -1, 1, -1, 1);
	free_vui8matrix(vX, -1, 1, -1, 1);
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
	max = (1 << (size * size)); // 2^25

	// Loop 1: Get permutations :
	printf("Implementation test : %s\n", dilation_set->func_name);
	for (perm = 0; perm < max; perm++) {
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
bool vec_morpho_produces_one(struct morpho_set *mset, vuint8** vW)
{
	vuint8 vX[5][3] = {{_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					   {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					   {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					   {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					   {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},};

	vuint8 vTempX[5][3] = {{_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					       {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					       {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					       {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},
					       {_mm_setzero_si128(),_mm_setzero_si128(),_mm_setzero_si128()},};	

	vuint8 *vY[5] = {vX[0] + 2, vX[1] + 2, vX[2] + 2, vX[3] + 2, vX[4] + 2};
	vuint8 *vTempY[5] = {vTempX[0] + 2, vTempX[1] + 2, vTempX[2] + 2, vTempX[3] + 2, vTempX[4] + 2};
	vuint8 **vZ = vY + 2;
	vuint8 **vTempZ = vTempY + 2;
	uint8 *z = (uint8*)&vZ[0][0];

	mset->vec_morpho_func(vW, 0, 0, 0, 0 /*, vTempZ*/, vZ);
	// int max = (1 << 9); // 2^9 
    
    // for (int i = 0; i < max; i++){
        
	// printf("vW:\n");
	// for (int row = -1; row <= 1; row++)
	// {
	// 	uint8 *p;
	// 	int i = 0;
	// 	for (int v = -1; v <= 1; v++)  {
	// 		p = (uint8*)&vW[row][v];
	// 		for(i=0; i<16; i++)
	// 			printf("%u", p[i]);
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");
	// printf("Z = %u\n", *z);
	// printf("vZ:\n");
	// for (int row = -1; row <= 1; row++)
	// {
	// 	uint8 *p;
	// 	int i = 0;
	// 	for (int v = -1; v <= 1; v++)  {
	// 		p = (uint8*)&vZ[row][v];
	// 		for(i=0; i<16; i++)
	// 			printf("%u", p[i]);
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");

    // printf("%s\n", z[0] == 1 ?"True" : "False");
	return z[0] == 1;
}
	// uint8 *x;
	// x = &X[1][1];
	// x[0]= W[0][0];
	// x[1]= W[0][1];
	// x[2]= W[0][2];
	// x[3]= W[0][3];
	// x[4]= W[0][4];



	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  0  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  10  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  210  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  3210  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  43210  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15   4321  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15    432  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15     43  */
	// /*0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15      4  */
	
// }



void test_intergration(uint8 **image, long nrl, long nrh, long ncl, long nch, const char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets,int nb_sets, bool logging)
{
      
    long temp_nrh, temp_nch;
    uint8 **X, **Y, **Z,**temp_buffer;

	assert((nrh - nrl + 1) > 20 && (nch - ncl + 1) > 20);
    // Test Input
    X = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// Middle Buffer
    temp_buffer = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	
    // Full zero intialization & copy image 
    memset_ui8matrix(X          , 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	memset_ui8matrix(temp_buffer, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
    copy_ui8matrix_ui8matrix(image, nrl, nrh, ncl, nch, X);
	
    for(temp_nrh = nrh - 20; temp_nrh < nrh + 1; temp_nrh++){
        for(temp_nch = nch - 20; temp_nch < nch + 1; temp_nch++){
			// Test Output
			Y = ui8matrix(nrl, temp_nrh, ncl, temp_nch);
			// Valid Output
			Z = ui8matrix(nrl, temp_nrh, ncl, temp_nch);
			for (int i = 0; i < nb_sets; i++) {
				if (morpho_sets[i].instr_type == SCALAR) {
				// if ((temp_nch + 1) % 10 == 0) {
					printf("Integration test : "LALIGNED_STR" (%-30s)\n", morpho_sets[i].func_name, filename);
					printf("\tTesting for %ld x %ld\n", temp_nch + 1, temp_nrh + 1);
				// }

					memset_ui8matrix(Y, 0, nrl, temp_nrh, ncl, temp_nch); 
					memset_ui8matrix(Z, 0, nrl, temp_nrh, ncl, temp_nch); 
					naive_morpho_set->morpho_func(X, nrl, temp_nrh, ncl, temp_nch, temp_buffer, Z);
					morpho_sets[i].morpho_func(   X, nrl, temp_nrh, ncl, temp_nch, temp_buffer, Y);

					if (logging) {
						printf("%ld %ld %ld %ld\n", nrl, temp_nrh, ncl, temp_nch);
						display_ui8matrix(X, nrl, temp_nrh, ncl, temp_nch, "%03u ", "Input");
						display_ui8matrix(Y, nrl, temp_nrh, ncl, temp_nch, "%03u ", morpho_sets[i].func_name);
						display_ui8matrix(Z, nrl, temp_nrh, ncl, temp_nch, "%03u ", naive_morpho_set->func_name);
					}
				// for (long row = nrl; row < temp_nrh; row++)
					assert(!memcmp_ui8matrix(Y, Z, nrl, temp_nrh, ncl, temp_nch));
				// if ((temp_nch + 1) % 10 == 0) 
					printf("\tTest passed\n");
				}
			}
			free_ui8matrix(Y, nrl, temp_nrh, ncl, temp_nch); 
			free_ui8matrix(Z, nrl, temp_nrh, ncl, temp_nch); 
		}
    }
	free_ui8matrix(temp_buffer, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
    free_ui8matrix(X, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
}
void test_packed_intergration(char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, const int nb_implementations, bool logging)
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
                if (logging) {
					printf("%ld %ld %ld %ld\n", nrl, temp_nrh, ncl, temp_nch);
                    display_ui8matrix(X, nrl, temp_nrh, ncl, temp_nch, "%u ", "Input");
                    display_ui8matrix(Y, nrl, temp_nrh, ncl, temp_nch, "%u ", morpho_sets[i].func_name);
                    display_ui8matrix(Z, nrl, temp_nrh, ncl, temp_nch, "%u ", naive_morpho_set->func_name);
                }
                assert(!memcmp_ui8matrix(Y, Z, nrl, temp_nrh, ncl, temp_nch));
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


void  ui8matrix_permutation (uint8** m, long nrl, long nrh, long ncl, long nch, uint32 perm)
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
	uint8 *x = &m[nrl][ncl];
	for (long i = 0; i < (nrh - nrl + 1) * (nch - ncl + 1); i++)
		x[i] = (perm >> i) & 0x1;
}

void  vui8matrix_permutation (vuint8** vM, long nrl, long nrh, long ncl, long nch, uint32 perm)
{
	/**
	 * Convert the permutation to a matrix format.
	 * 
	 * For example: 
	 * perm = 0b111000111
	 * =>
	 * X = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},{1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}
	 *     {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}
	 *     {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},{1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}
	 **/
	int card = card_vuint8();
	assert(-card < ncl && nch < card * 2);
	
	int ncl0 = (card - 1) + (ncl < 0 ? ncl : 0) + 1;
	int nch0 = (card - 1);

	int ncl1 = (ncl >= 0    ? ncl : 0       );
	int nch1 = (nch <  card ? nch : card - 1);

	int ncl2 = (ncl >= card ? ncl - card : 0 );
	int nch2 = (nch >= card ? nch - card : 0 ) - 1;
	int row, col, row_cnt = 0, col_cnt = 0, ncol = (nch - ncl + 1);
	zero_vui8matrix(vM, nrl, nrh, -1, 1);
	uint8 *v;

	row_cnt = 0;
	for (row = nrl; row < nrh + 1; row++, row_cnt += ncol) {
		col_cnt = 0;
		v = (uint8 *)&vM[row][-1];

		for (col = ncl0; col < nch0 + 1; col++) 
			v[col] = (perm >> (row_cnt +col_cnt++)) & 0x1;

		v = (uint8 *)&vM[row][ 0];
		for (col = ncl1; col < nch1 + 1; col++) 
			v[col] = (perm >> (row_cnt +col_cnt++)) & 0x1;


		v = (uint8 *)&vM[row][ 1];
		for (col = ncl2; col < nch2 + 1; col++) 
			v[col] = (perm >> (row_cnt +col_cnt++)) & 0x1;
	}

	// int max = (1 << 9); // 2^9 
    
    // for (int i = 0; i < max; i++){
    //     printf("%x\n",i);
    //     vuint8 **m = vui8matrix(-1,1,-1,1);
    //     vui8matrix_permutation(m, -1, 1, -1, 1, i);
        

    //     for (int row = -1; row <= 1; row++)
    //     {
    //         uint8 *p;
    //         int i = 0;
    //         for (int v = -1; v <= 1; v++)  {
    //             p = (uint8*)&m[row][v];
    //             for(i=0; i<16; i++)
    //                 printf("%u", p[i]);
    //         }
    //         printf("\n");
    //     }
    //         printf("\n");
        
    //     getchar();
        
    //     free_vui8matrix(m, -1, 1, -1, 1);
    // }
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
					    (row + y) <= nrh && (col + x) <= nch) {
						Y[row + y][col + x] = xor_mask;
						Y[row + y][col + x] ^= xor_mask;
					}
			}
		}
	}
	// display_ui8matrix(Y, nrl, nrh, ncl, nch, "%u", "CHekcer");
	// getchar();
	return Y;
}

