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
uint8 **erosion(p_image img, p_struct_elem s) 
{
	// bord = border
	// th = top horizontal / bv = bottom horizontal
	// lv = left vertical / rh = right vertical
    long nthbord, nlvbord, nbhbord, nrvbord;
	uint8 **input, **output;
	long col, row, y, x;
	uint8 pixel;

	// Compute the size of borders (preduplication).
	nthbord = s->oriy;
	nbhbord = (s->nrow - 1) - s->oriy;
	nlvbord = s->orix;
	nrvbord = (s->ncol - 1) - s->orix;
	
	// Create a new matrix with edges 
	input = ui8matrix(img->nrl - nthbord, img->nrh + nbhbord, img->ncl - nlvbord, img->nch + nrvbord);
	// Copy the original image to the input matrix.
	copy_ui8matrix_ui8matrix(img->I, img->nrl, img->nrh, img->ncl, img->nch, input);	
	// Erode
	for (row = img->nrl; row < img->nrh + 1; row++){
		for (col = img->ncl; col < img->nch + 1; col++) {
			pixel = input[row][col];
			for (y = -nthbord; y < nbhbord + 1; y++ ) {
				for (x = -nlvbord; x < nrvbord + 1; x++ ) {
					pixel &= input[row + y][col + x];
				}
			}
			img->Omega[row][col] = pixel;
		}
	}
	// Save the image (debug)
	SavePGM_ui8matrix(img->Omega, img->nrl, img->nrh, img->ncl, img->nch, "output_erosion.pgm");
	// Free the input matrix.
	free_ui8matrix(input, img->nrl - nthbord, img->nrh + nbhbord, img->ncl - nlvbord, img->nch + nrvbord);
	return output;
}
uint8 **dilation(p_image img, p_struct_elem s) 
{
	// bord = border
	// th = top horizontal / bv = bottom horizontal
	// lv = left vertical / rh = right vertical
    long nthbord, nlvbord, nbhbord, nrvbord;
	uint8 **input, **output;
	long col, row, y, x;
	uint8 pixel;

	// Compute the size of borders (preduplication).
	nthbord = s->oriy;
	nbhbord = (s->nrow - 1) - s->oriy;
	nlvbord = s->orix;
	nrvbord = (s->ncol - 1) - s->orix;
	
	// Create a new matrix with edges 
	input = ui8matrix(img->nrl - nthbord, img->nrh + nbhbord, img->ncl - nlvbord, img->nch + nrvbord);
	// Copy the original image to the input matrix.
	copy_ui8matrix_ui8matrix(img->I, img->nrl, img->nrh, img->ncl, img->nch, input);	
	// Erode
	for (row = img->nrl; row < img->nrh + 1; row++){
		for (col = img->ncl; col < img->nch + 1; col++) {
			pixel = input[row][col];
			for (y = -nthbord; y < nbhbord + 1; y++ ) {
				for (x = -nlvbord; x < nrvbord + 1; x++ ) {
					pixel |= input[row + y][col + x];
				}
			}
			img->Omega[row][col] = pixel;
		}
	}
	// Save the image (debug)
	SavePGM_ui8matrix(img->Omega, img->nrl, img->nrh, img->ncl, img->nch, "output_dilation.pgm");
	// Free the input matrix.
	free_ui8matrix(input, img->nrl - nthbord, img->nrh + nbhbord, img->ncl - nlvbord, img->nch + nrvbord);
	return output;
}
#define TEST_IMAGE_CX 21
#define TEST_IMAGE_CY 17
void test_morpho()
{
	p_struct_elem s = create_structuring_element(1,1,3,3);
	// for (int i = 0; i < s->nrow; i++) {x`
	// 	for (int j = 0; j < s->ncol; j++) {
	// 		s->m[i][j] = mask3x3_plus[i][j];
	// 	}
	// }
	// display_ui8matrix(s->m, 0, s->nrow - 1, 0, s->ncol - 1, "%u ", "mask3x3_plus");
	erosion(create_image("./morpho_test.pgm"), s);
	dilation(create_image("./morpho_test.pgm"), s);
}