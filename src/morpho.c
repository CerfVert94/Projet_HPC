#include "nrdef.h"
#include "nrutil.h"
#include <img.h>
#include <stdlib.h>
#include <stdio.h>
//#include <malloc.h>
#include <morpho.h>


	// uint8 mask3x3_plus[3][3] = {{0,1,0},
	// 							{1,1,1},
	// 							{0,1,0}};


	// uint8 mask3x3_ltri[3][3] = {{1,0,0},
	// 							{1,1,0},
	// 							{1,1,1}};

	// uint8 mask3x3_htri[3][3] = {{1,1,1},
	// 							{0,1,1},
	// 							{0,0,0}};

	// uint8 mask3x3_flat[3][3] = {{1,1,1},
	// 							{1,1,1},
	// 							{1,1,1}};

	// uint8 mask2x2_ltri[2][2] = {{1,0},
	// 							{1,1}};							

	// uint8 mask2x2_htri[2][2] = {{1,1},
	// 							{0,1}};							

	// uint8 mask2x2_flat[2][2] = {{1,1},
	// 							{0,1}};							
		

p_struct_elem create_structuring_element(long orix, long oriy, long nrow, long ncol)
{
	p_struct_elem s;

	s = (p_struct_elem) malloc(sizeof(struct_elem));

	// Create a nrow x ncol uint8 matrix.
	// s->m = ui8matrix(0, nrow - 1, 0, ncol - 1);

	// Define dimension.
	s->nrow = nrow;
	s->ncol = ncol;
	// Define origin.
	s->orix= orix;
	s->oriy= oriy;
}

// Without optimisation
uint8** ui8matrix_dilation(uint8** input, long nrl, long nrh, long ncl, long nch, p_struct_elem s)
{
	// bord = border
	// th = top horizontal / bv = bottom horizontal
	// lv = left vertical / rh = right vertical
    long nthbord, nlvbord, nbhbord, nrvbord;
	uint8 **bord_input, **output;
	long col, row, y, x;
	uint8 pixel;

	// Compute the size of borders (preduplication).
	nthbord = s->oriy;
	nbhbord = (s->nrow - 1) - s->oriy;
	nlvbord = s->orix;
	nrvbord = (s->ncol - 1) - s->orix;
	
	// Create a new matrix with edges / output matrix
	bord_input  = ui8matrix(nrl - nthbord, nrh + nbhbord, ncl - nlvbord, nch + nrvbord);
	output 		= ui8matrix(nrl, nrh, ncl, nch);
	
	// Copy the original image to the input matrix.
	copy_ui8matrix_ui8matrix(input, nrl, nrh, ncl, nch, bord_input);
	
	// Erode
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			pixel = bord_input[row][col];
			for (y = -nthbord; y < nbhbord + 1; y++ ) {
				for (x = -nlvbord; x < nrvbord + 1; x++ ) {
					pixel |= bord_input[row + y][col + x];
				}
			}
			output[row][col] = pixel;
		}
	}
	
	// Free the bordered matrix.
	free_ui8matrix(bord_input, nrl - nthbord, nrh + nbhbord, ncl - nlvbord, nch + nrvbord);
	return output;
}
uint8** ui8matrix_erosion(uint8** input, long nrl, long nrh, long ncl, long nch, p_struct_elem s)
{
	// bord = border
	// th = top horizontal / bv = bottom horizontal
	// lv = left vertical / rh = right vertical
    long nthbord, nlvbord, nbhbord, nrvbord;
	uint8 **bord_input, **output;
	long col, row, y, x;
	uint8 pixel;

	// Compute the size of borders (preduplication).
	nthbord = s->oriy;
	nbhbord = (s->nrow - 1) - s->oriy;
	nlvbord = s->orix;
	nrvbord = (s->ncol - 1) - s->orix;
	printf("%ld %ld %ld %ld\n", nrl, nrh, ncl, nch);
	// Create a new matrix with edges / output matrix
	bord_input  = ui8matrix(nrl - nthbord, nrh + nbhbord, ncl - nlvbord, nch + nrvbord);
	output 		= ui8matrix(nrl, nrh, ncl, nch);

	// Copy the original image to the input matrix.
	copy_ui8matrix_ui8matrix(input, nrl, nrh, ncl, nch, bord_input);	
	// Dilate
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			pixel = bord_input[row][col];
			for (y = -nthbord; y < nbhbord + 1; y++ ) {
				for (x = -nlvbord; x < nrvbord + 1; x++ ) {
					pixel &= bord_input[row + y][col + x];
				}
			}
			output[row][col] = pixel;
		}
	}
	// Free the bordered matrix.
	free_ui8matrix(bord_input, nrl - nthbord, nrh + nbhbord, ncl - nlvbord, nch + nrvbord);
	return output;
}


void image_dilation(p_image img, p_struct_elem s)
{
	uint8 **output;
	long nrl, nrh, ncl, nch;
	
	nrl = img->nrl; nrh = img->nrh;
	ncl = img->ncl;	nch = img->nch;

	output = ui8matrix_dilation(img->I, nrl, nrh, ncl, nch, s);	
	copy_ui8matrix_ui8matrix(output, nrl, nrh, ncl, nch, img->Omega);
	// Save the image (debug)
	SavePGM_ui8matrix(img->Omega, nrl, nrh, ncl, nch,"output_dilation.pgm");
	free_ui8matrix(output, nrl, nrh, ncl, nch);
}
void image_erosion(p_image img, p_struct_elem s)
{
	uint8 **output;
	long nrl, nrh, ncl, nch;
	
	nrl = img->nrl; nrh = img->nrh;
	ncl = img->ncl;	nch = img->nch;

	output = ui8matrix_erosion(img->I, nrl, nrh, ncl, nch, s);
	copy_ui8matrix_ui8matrix(output, nrl, nrh, ncl, nch, img->Omega);
	// Save the image (debug)
	SavePGM_ui8matrix(img->Omega, nrl, nrh, ncl, nch,"output_erosion.pgm");
	free_ui8matrix(output, nrl, nrh, ncl, nch);
}
void test_morpho()
{
	// const int TEST_CX = 3;
	// const int TEST_CY = 3;
	// long y, x;
	p_struct_elem s = create_structuring_element(1,1,3,3);
	// const uint8 input3x3[3][3]  = {{1,1,1},
	// 							   {1,1,1},
	// 							   {1,1,1}};
	// const uint8 output3x3[3][3] = {{0,0,0},
	// 							   {0,1,0},
	// 							   {0,0,0}};
	// uint8 **res;
	// res = erosion(input3x3, s);

	image_erosion(create_image("morpho_test.pgm"), s);
	image_dilation(create_image("morpho_test.pgm"), s);
}